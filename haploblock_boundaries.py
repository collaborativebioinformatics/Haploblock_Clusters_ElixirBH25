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
    logger.info("Found %i positions with recombination rates > 10 x average", len(high_rates))
    
    # add first haploblock
    haploblock_boundaries.append((1, high_rates_positions[0]))
    
    for i in range(1, len(high_rates)):
        start = high_rates_positions[i - 1]
        end = high_rates_positions[i]
        haploblock_boundaries.append((start, end))
    
    # add last haploblock
    haploblock_boundaries.append((high_rates_positions[-1], positions[-1]))
    
    return(haploblock_boundaries)


def haploblocks_to_TSV(haploblock_boundaries):
    '''
    Print haploblock boundaries to stdout in TSV format, 2 columns: start end

    arguments:
    - haploblock_boundaries: list of tuples with haploblock boundaries (start, end)
    '''
    # header
    print("START\tEND")

    for (start, end) in haploblock_boundaries:
        print(str(start) + "\t" + str(end))


def main(recombination_file, chromosome):

    logger.info("Parsing recombination file")
    haploblock_boundaries = parse_recombination_rates(recombination_file, chromosome)
    
    logger.info("Printing haploblock boundaries")
    haploblocks_to_TSV(haploblock_boundaries)


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

    args = parser.parse_args()

    try:
        main(recombination_file=args.recombination_file,
             chromosome=args.chr)

    except Exception as e:
        # details on the issue should be in the exception name, print it to stderr and die
        sys.stderr.write("ERROR in " + script_name + " : " + repr(e) + "\n")
        sys.exit(1)