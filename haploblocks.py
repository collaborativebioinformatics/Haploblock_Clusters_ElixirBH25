import os
import sys
import logging
import argparse
import pathlib


# set up logger, using inherited config, in case we get called as a module
logger = logging.getLogger(__name__)


def parse_recombination_rates(recombination_file, chromosome):
    """
    Parses recombination rates from Halldorsson et al., 2019

    arguments:
    - recombination_file with 7-line header
      and 5 columns: Chr Begin End cMperMb cM

    returns:
    - haploblock_boundaries: list of tuples with haploblock boundaries (start, end)
    """
    if not chromosome.startswith("chr"):
        chromosome = "chr" + str(chromosome)
    
    haploblock_boundaries = []  # start at pos==1
    positions = []
    rates = []
    high_rates = []
    high_rates_positions = []

    try:
        f = open(recombination_file, 'r')
    except Exception as e:
        logger.error("Opening provided recombination file %s: %s", recombination_file, e)
        raise Exception("Cannot open provided recombination file")

    for line in f:
        # skip header
        if line.startswith("#"):
            continue

        split_line = line.rstrip().split('\t')

        if len(split_line) != 5:
            logger.error("Recombination file %s has bad line (not 5 tab-separated fields): %s",
                         recombination_file, line)
            raise Exception("Bad line in the recombination file")

        (chr, start, end, rate, cM) = split_line

        if chr == chromosome:
            positions.append(int(start))
            rates.append(float(rate))
    
    assert len(positions) == len(rates)

    # find positions with recombination rate > 10*avg
    avg_rate = sum(rates) / len(rates)

    for i in range(len(rates)):
        if rates[i] > 10 * avg_rate:
            high_rates.append(rates[i])
            high_rates_positions.append(positions[i])
    
    # add first haploblock
    haploblock_boundaries.append((1, high_rates_positions[0]))
    
    for i in range(1, len(high_rates)):
        start = high_rates_positions[i - 1]
        end = high_rates_positions[i]
        haploblock_boundaries.append((start, end))
    
    # add last haploblock
    haploblock_boundaries.append((high_rates_positions[-1], positions[-1]))

    logger.info("Found %i haploblocks", len(haploblock_boundaries))
    
    return(haploblock_boundaries)


def generate_haploblock_hashes(haploblock_boundaries):
    """
    Generate a unique identifier, "haploblock hash", for every haploblock,
    ie an integer with len(haploblock_boundaries) digits

    returns:
    - haploblock2hash: dict, key=(start, end), value=haploblock hash
    """
    haploblock2hash = {}
    
    idx = len(haploblock_boundaries) - 1  # start hashes from the end
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
    haploblock_boundaries = parse_recombination_rates(recombination_file, chr)

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