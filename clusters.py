import os
import sys
import logging
import argparse
import pathlib
import subprocess

import data_parser


# set up logger, using inherited config, in case we get called as a module
logger = logging.getLogger(__name__)


def calculate_mmseq_params(variant_counts_file):
    """
    Parses variant counts file with 4 columns: start, end, mean, stdev
    mean, stdev of the number of variants per haploblock
    calculates min identity

    arguments:
    - variant_counts_file
    returns:
    - haploblock2min_id: dict, key=(start, end), value=min identity
    - haploblock2cov_fraction: dict, key=(start, end), value=coverage fraction
    """
    haploblock2min_id = {}
    haploblock2cov_fraction = {}

    try:
        f = open(variant_counts_file, 'r')
    except Exception as e:
        logger.error("Opening provided variant counts file %s: %s", variant_counts_file, e)
        raise Exception("Cannot open provided variant counts file")
    
    # skip header
    line = f.readline()
    if not line.startswith("START\t"):
        logging.error("variant counts file %s is headerless? expecting headers but got %s",
                      variant_counts_file, line)
        raise Exception("variant counts file problem")

    for line in f:
        split_line = line.rstrip().split('\t')

        if len(split_line) != 4:
            logger.error("variant counts file %s has bad line (not 4 tab-separated fields): %s",
                         variant_counts_file, line)
            raise Exception("Bad line in the variant counts file")
        
        (start, end, mean, stdev) = split_line

        haploblock_len = int(end) - int(start)
        haploblock2min_id[(start, end)] = 1 - (float(mean) / haploblock_len)
        haploblock2cov_fraction[(start, end)] = 1 - (682 / haploblock_len)

    return(haploblock2min_id, haploblock2cov_fraction)


def compute_clusters(input_fasta, out, min_seq_id, cov_fraction, cov_mode,
                     chr, start, end):
    """
    Run MMSeqs2

    arguments:
    - input_fasta: merged phased fasta file for haploblock
    """

    output_prefix = os.path.join(out, "clusters", f"chr{chr}_{start}-{end}")

    subprocess.run(["mmseqs", "easy-cluster",
                    input_fasta,
                    output_prefix,
                    os.path.join(out, "tmp"),
                    "--min-seq-id", str(min_seq_id),
                    "-c", str(cov_fraction),
                    "--cov-mode", str(cov_mode),
                    "--remove-tmp-files", "1"],
                    check=True)


def main(boundaries_file, merged_consensus_dir, variant_counts_file, chr, out, cov_mode):

    if os.path.exists(os.path.join(out, "clusters")):
        logger.error(f"Output directory {os.path.join(out)} exists, please remove it")
        raise Exception("Output directory exists")
    os.mkdir(os.path.join(out, "clusters"))


    logger.info("Parsing haploblock boundaries")
    haploblock_boundaries = data_parser.parse_haploblock_boundaries(boundaries_file)
    logger.info("Found %i haploblocks", len(haploblock_boundaries))

    logger.info("Parsing variant counts file")
    (haploblock2min_id, haploblock2cov_fraction) = calculate_mmseq_params(variant_counts_file)

    logger.info("Calculating clusters")
    for (start, end) in haploblock_boundaries:
        input_fasta = os.path.join(merged_consensus_dir, f"chr{chr}_region_{start}-{end}.fa")
        min_seq_id = haploblock2min_id[(start, end)]
        cov_fraction = haploblock2cov_fraction[(start, end)]

        compute_clusters(input_fasta, out, min_seq_id, cov_fraction, cov_mode,
                         chr, start, end)


if __name__ == "__main__":
    script_name = os.path.basename(sys.argv[0])
    # configure logging, sub-modules will inherit this config
    logging.basicConfig(format='%(asctime)s %(levelname)s %(name)s: %(message)s',
                        datefmt='%Y-%m-%d %H:%M:%S',
                        level=logging.DEBUG)
    # set up logger: we want script name rather than 'root'
    logger = logging.getLogger(script_name)

    parser = argparse.ArgumentParser(
        prog=script_name,
        description="TODO"
    )
    
    parser.add_argument('--boundaries_file',
                        help='Path to boundaries file generated from Halldorsson et al., 2019',
                        type=pathlib.Path,
                        required=True)
    parser.add_argument('--merged_consensus_dir',
                        help='Path to directory with consensus haploblock phased sequences fasta files per haploblock',
                        type=pathlib.Path,
                        required=True)
    parser.add_argument('--variant_counts',
                        help='Path to file with the mean and average of the number of variants per haploblock',
                        type=pathlib.Path,
                        required=True)
    parser.add_argument('--chr',
                        help='chromosome',
                        type=str,
                        required=True)
    parser.add_argument('--out',
                        help='Path to output directory for clusters',
                        type=pathlib.Path,
                        required=True)
    parser.add_argument('--cov_mode',
                        help='coverage mode (optional)',
                        type=int,
                        required=False,
                        default="0")

    args = parser.parse_args()

    try:
        main(boundaries_file=args.boundaries_file,
             merged_consensus_dir=args.merged_consensus_dir,
             variant_counts_file=args.variant_counts,
             chr=args.chr,
             out=args.out,
             cov_mode=args.cov_mode)

    except Exception as e:
        # details on the issue should be in the exception name, print it to stderr and die
        sys.stderr.write("ERROR in " + script_name + " : " + repr(e) + "\n")
        sys.exit(1)