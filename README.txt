**Simple VCF parser**

Not terribly efficient, as it just reads through lines in the VCF and tallies genotype counts.
But does not require loading a large genotype matrix and has some features for filtering.

Usage:

    python parse_vcf.py --input input.vcf(.gz) --output allele_counts.bp --positions positions.txt/bp --ancestral_sequence anc.fa(.gz) --verbose 1000 --sample_file samples.txt

See argparse in script for more information about usage.

