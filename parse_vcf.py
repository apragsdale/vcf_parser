# Script to take input options and run through a VCF to compute allele counts for SFS.
# Input options are samples file, positions file, CSQ (?), multiallelic, SNPs only,
# ancestral sequence for polarization.
# Input file must end in .vcf, .vcf.gz, or .vcf.bgz, and automatically detects
# if file is gzip compressed or not.
#


import sys
import gzip
import argparse
from collections import defaultdict
from Bio import SeqIO
import pickle


def make_parser():
    ADHF = argparse.ArgumentDefaultsHelpFormatter
    parser = argparse.ArgumentParser("get_sfs.py", formatter_class=ADHF)
    parser.add_argument(
        "--input", "-i", type=str, help="Input VCF file (.vcf, .vcf.gz, .vcf.bgz)."
    )
    parser.add_argument(
        "--output",
        "-o",
        type=str,
        help="Output file (pickled dictionary of allele counts).",
    )
    optional = parser.add_argument_group("Optional")
    optional.add_argument(
        "--sample_file",
        "-s",
        type=str,
        default=None,
        help="File assigning samples to populations.",
    )
    optional.add_argument(
        "--positions",
        type=str,
        default=None,
        help="File of positions: text file with lines 'chr pos', or pickled dictionary.",
    )
    optional.add_argument(
        "--keep_multiallelic",
        action="store_true",
        help="If keep_multiallelic is set, we sum all derived alleles together.",
    )
    optional.add_argument(
        "--all_variants",
        action="store_true",
        help="If all_variants is set, we consider all allele types. Otherwise, only SNPs.",
    )
    optional.add_argument(
        "--ancestral_sequence",
        "-a",
        type=str,
        default=None,
        help="Ancestral sequence in fasta format (.fasta, .fa, .fasta.gz, .fa.gz).",
    )
    optional.add_argument(
        "--filter", "-f", type=int, default=None, help="Minimum genome quality to keep"
    )
    optional.add_argument(
        "--verbose",
        "-v",
        type=int,
        default=None,
        help="Spacing for reporting progress along chromosome",
    )
    return parser


def get_pop_cols(header_line, sample_file):
    """
    header_line: the line that starts with #CHROM, decoded.
    sample_file: each line is "sample{white-space}pop"
    """
    if sample_file is None:
        return {"ALL": list(range(9, len(header_line.split())))}
    else:
        pops = defaultdict(list)
        for line in open(sample_file, "r"):
            pops[line.split()[1]].append(line.split()[0])
        pop_cols = defaultdict(list)
        for pop in pops.keys():
            for sid in pops[pop]:
                try:
                    pop_cols[pop].append(header_line.split().index(sid))
                except ValueError:
                    raise ValueError(
                        "{sid} in population {pop} is not in header samples"
                    )
        return pop_cols


def get_pop_genotypes(line_split, pop_cols):
    gts = [line_split[i].split(":")[0] for i in pop_cols]
    return gts


def count_genotypes(gt, multiallelic=False, anc_idx=0):
    """
    gt: list of genotypes, stripped of 
    Given list of ./., returns allele count.
    If multiallelic is True, we have to sum all derived variants.
    """
    if multiallelic:
        raise ValueError("Multiallelic not implemented")
    else:
        ac = (
            gt.count("0/1")
            + gt.count("0|1")
            + gt.count("1/0")
            + gt.count("1|0")
            + 2 * gt.count("1/1")
            + 2 * gt.count("1|1")
        )
        n_hap = 2 * (
            gt.count("0/0")
            + gt.count("0|0")
            + gt.count("0/1")
            + gt.count("0|1")
            + gt.count("1/0")
            + gt.count("1|0")
            + gt.count("1/1")
            + gt.count("1|1")
        )
        if anc_idx == 1:
            ac = n_hap - ac
    return ac, n_hap


if __name__ == "__main__":
    parser = make_parser()
    args = parser.parse_args(sys.argv[1:])

    if args.verbose is not None:
        print("Processing:")
        print("input:", args.input)
        print("output:", args.output)
    
    fin = args.input
    fout = args.output

    if fin.endswith("vcf"):
        open_func = open
    elif fin.endswith("vcf.gz") or fin.endswith("vcf.bgz"):
        open_func = gzip.open
    else:
        raise ValueError("Invalid input file - must end with vcf or vcf.gz/bgz")

    if args.positions is not None:
        # get dict for positions
        if args.positions.endswith(".txt"):
            pos_dict = defaultdict(list)
            for line in open(args.positions):
                pos_dict[int(line.split()[0])].append(int(line.split()[1]))
        else:
            pos_dict = pickle.load(open(args.positions, "rb"))

    if args.ancestral_sequence is not None:
        # load ancestral sequence fasta
        if args.ancestral_sequence.endswith(".gz"):
            anc_f = gzip.open(args.ancestral_sequence, "rt")
        elif args.ancestral_seq.endswith(".fa") or args.ancestral_sequence.endswith(
            ".fasta"
        ):
            anc_f = open(args.ancestral_sequence, "rt")
        else:
            raise ValueError("not fasta or gzipped fasta?")

        for record in SeqIO.parse(anc_f, "fasta"):
            anc_seq = record.seq
        anc_f.close()

    snps = ["A", "C", "G", "T"]

    # keep data on kept and skipped variants
    pos_kept = 0
    no_match = 0
    match_ref = 0
    match_alt = 0
    skipped_multi = 0

    allele_counts = defaultdict(lambda: defaultdict(int))

    with open_func(fin, "rb") as f:
        for line in f:
            l = line.decode()
            if l.startswith("#"):
                if l.startswith("#CHROM"):
                    pop_cols = get_pop_cols(l, args.sample_file)
                    pops = [p for p in pop_cols.keys()]
            else:
                lsplit = l.split()
                chrom, pos, _, ref, alt, filt = lsplit[:6]
                chrom = int(chrom)
                pos = int(pos)

                if args.positions is not None:
                    if pos not in pos_dict[chrom]:
                        continue

                if args.keep_multiallelic is False:
                    if len(alt.split(",")) > 1:
                        skipped_multi += 1
                        continue

                if args.all_variants is False:
                    if ref not in snps:
                        continue
                    for a in alt.split(","):
                        if a not in snps:
                            continue

                if args.ancestral_sequence is not None:
                    anc_idx = -1
                    if anc_seq[pos - 1].upper() == ref:
                        anc_idx = 0
                        match_ref += 1
                    else:
                        for ii, a in enumerate(alt.split(",")):
                            if anc_seq[pos - 1].upper() == a:
                                anc_idx = ii + 1
                                match_alt += 1
                    if anc_idx == -1:
                        no_match += 1
                        continue
                else:
                    match_ref = 1
                    anc_idx = 0

                acs_ns = [
                    count_genotypes(
                        get_pop_genotypes(lsplit, pop_cols[p]),
                        multiallelic=args.keep_multiallelic,
                        anc_idx=anc_idx,
                    )
                    for p in pops
                ]

                acs = tuple([_[0] for _ in acs_ns])
                ns = tuple([_[1] for _ in acs_ns])

                allele_counts[ns][acs] += 1

                pos_kept += 1
                if args.verbose is not None:
                    if pos_kept % args.verbose == 0:
                        print(f"kept {pos_kept} sites, at position {pos}")

    data = {}
    data["allele_counts"] = dict(allele_counts)
    data["stats"] = {
        "kept": pos_kept,
        "no_match": no_match,
        "match_ref": match_ref,
        "match_alt": match_alt,
    }
    
    if args.verbose is not None:
        print("final stats:", data["stats"])
    with open(args.output, "wb+") as fout:
        pickle.dump(data, fout)
