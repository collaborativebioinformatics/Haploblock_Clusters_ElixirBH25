import os
import sys
import logging
import argparse
import pathlib

import data_parser

# set up logger, using inherited config, in case we get called as a module
logger = logging.getLogger(__name__)


def generate_haploblock_hashes(haploblock_boundaries):
    """
    Generate a unique identifier, "haploblock hash", for every haploblock,
    ie an integer with len(haploblock_boundaries) digits

    returns:
    - haploblock2hash: dict, key=(start, end), value=haploblock hash
    """
    haploblock2hash = {}
    
    idx = len(haploblock_boundaries) - 1  # populating cluster hashes from the end
    for (start, end) in haploblock_boundaries:
        hashList = ["0"] * len(haploblock_boundaries)
        hashList[idx] = "1"
        hash = "".join(hashList)
        haploblock2hash[(start, end)] = hash
        idx -= 1
    
    return(haploblock2hash)


def haploblocks_to_TSV(haploblock_boundaries, chr, out):
    '''
    Save haploblock boundaries to a TSV file, 2 columns: start end

    arguments:
    - haploblock_boundaries: list of tuples with haploblock boundaries (start, end)
    '''
    with open(os.path.join(out, f"haploblock_boundaries_chr{chr}.tsv"), 'w') as f:
        # header
        f.write("START\tEND\n")
        for (start, end) in haploblock_boundaries:
            f.write(str(start) + "\t" + str(end) + "\n")


def haploblock_hashes_to_TSV(haploblock2hash, chr, out):
    '''
    Save haploblock hashes to a TSV file, " columns: start end hash

    arguments:
    - haploblock2hash: dict, key=(start, end), value=haploblock hash
    '''
    with open(os.path.join(out, f"haploblock_hashes_chr{chr}.tsv"), 'w') as f:
        # header
        f.write("START\tEND\tHASH\n")
        for (start, end) in haploblock2hash:
            f.write(str(start) + "\t" + str(end) +  "\t" + haploblock2hash[(start, end)] + "\n")


def main(recombination_file, chr, out):

    logger.info("Parsing recombination file")
    haploblock_boundaries = data_parser.parse_recombination_rates(recombination_file, chr)

    logger.info("Generating haploblock hashes")
    haploblock2hash = generate_haploblock_hashes(haploblock_boundaries)
    
    logger.info("Saving haploblock boundaries")
    haploblocks_to_TSV(haploblock_boundaries, chr, out)

    logger.info("Saving haploblock hashes")
    haploblock_hashes_to_TSV(haploblock2hash, chr, out)


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
    
    parser.add_argument('--recombination_file',
                        help='Path to recombination file from Halldorsson et al., 2019',
                        type=pathlib.Path,
                        required=True)
    parser.add_argument('--chr',
                        help='chromosome',
                        type=str,
                        required=True)
    parser.add_argument('--out',
                        help='Path to output folder',
                        type=pathlib.Path,
                        required=True)

    args = parser.parse_args()

    try:
        main(recombination_file=args.recombination_file,
             chr=args.chr,
             out=args.out)

    except Exception as e:
        # details on the issue should be in the exception name, print it to stderr and die
        sys.stderr.write("ERROR in " + script_name + " : " + repr(e) + "\n")
        sys.exit(1)