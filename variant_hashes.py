import os
import sys
import logging
import argparse
import pathlib

import data_parser

# set up logger, using inherited config, in case we get called as a module
logger = logging.getLogger(__name__)


def generate_cluster_hashes(clusters):
    """
    Generate unique cluster hashes for each cluster,
    ie strings of len(clusters) with 1/0s
    
    arguments:
    - clusters: list of cluster IDs

    returns:
    - cluster2hash: dict, key=clusterID, value=hash
    """
    cluster2hash = {}

    idx = len(clusters) - 1  # populating cluster hashes from the end
    for cluster in clusters:
        hashList = ["0"] * len(clusters)
        hashList[idx] = "1"
        hash = "".join(hashList)
        cluster2hash[cluster] = hash
        idx -= 1

    return(cluster2hash)


def generate_individual_hashes(individual2cluster, cluster2hash, variant2hash, haploblock2hash, chr_hash):
    """
    Generate individual hashes, ie 64-character strings of 0/1s, each contains:
        strand hash: 4 char
        chromosome hash: 10 chars
        haploblock hash: len(haploblocks) chars
        cluster hash: len(clusters) chars
        variant hash: len(variants of interest) chars

    arguments:
    - individual2cluster: dict, key=individual, value=unique clusterID
    - cluster2hash: dict, key=clusterID, value=hash
    - variant2hash: dict, key=individual, key=hash
    - haploblock2hash: dict, key=(start, end), key=hash
    - chr_hash: 10-digit

    returns:
    - individual2hash: dict, key=individual, value=hash
    """
    individual2hash = {}

    for individual in individual2cluster:
        # strand hash
        # individual is in format: NA18524_chr6_region_31480875-31598421_hap0
        # this should be fixed in the future, for now taking the last char from individual as strand (can be 0 or 1)
        strand = individual[-1]
        strand_hash = ""
        if strand == "0":
            strand_hash = "0001"
        else:
            strand_hash = "0010"

        # haploblock_hash
        # individual is in format: NA18524_chr6_region_31480875-31598421_hap0
        # this should be fixed in the future, for now taking start-end
        individual_split = individual.split("_")
        (start, end) = individual_split[3].split("-")
        haploblock_hash = haploblock2hash[(start, end)]

        # cluster hash
        cluster = individual2cluster[individual]
        cluster_hash = cluster2hash[cluster]

        # variant hash
        variant_hash = variant2hash[individual]

        # individual hash
        hash = strand_hash + chr_hash + haploblock_hash + cluster_hash + variant_hash
        individual2hash[individual] = hash

    return(individual2hash)


def hashes_to_TSV(individual2hash, out):
    '''
    Save individual to a TSV file, 2 columns: INDIVIDUAL HASH

    arguments:
    - individual2hash: dict, key=individual, value=hash
    '''
    with open(os.path.join(out, f"individual_hashes.tsv"), 'w') as f:
        # header
        f.write("INDIVIDUAL\tHASH\n")
        for individual in individual2hash:
            f.write(individual + "\t" + individual2hash[individual] + "\n")


def main(clusters_file, variant_hashes, haploblock_hashes, out):

    logger.info("Parsing clusters")
    # this should be fixed, beacuse we will have one clusters file per haploblock, not one clusters file in total
    (individual2cluster, clusters) = data_parser.parse_clusters(clusters_file)
    logger.info("Found %i clusters with %i individuals in total", len(clusters), len(individual2cluster))

    logger.info("Parsing variant hashes")
    variant2hash = data_parser.parse_variant_hashes(variant_hashes)
    
    logger.info("Parsing haploblock hashes")
    haploblock2hash = data_parser.parse_haploblock_hashes(haploblock_hashes)

    logger.info("Generating individual hashes")
    # for now hardcode chromosome hash
    chr_hash = "0000100000"  # chr6  
    cluster2hash = generate_cluster_hashes(clusters)
    individual2hash = generate_individual_hashes(individual2cluster, cluster2hash, variant2hash, haploblock2hash, chr_hash)

    logger.info("Saving individual hashes")
    hashes_to_TSV(individual2hash, out)


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
    
    parser.add_argument('--clusters',
                        help='Path to clusters file generated with MMSeqs2',
                        type=pathlib.Path,
                        required=True)
    parser.add_argument('--variant_hashes',
                        help='Path to file with variant hashes',
                        type=pathlib.Path,
                        required=True)
    parser.add_argument('--haploblock_hashes',
                        help='Path to file with haploblock hashes',
                        type=pathlib.Path,
                        required=True)
    parser.add_argument('--out',
                        help='Path to output folder',
                        type=pathlib.Path,
                        required=True)

    args = parser.parse_args()

    try:
        main(clusters_file=args.clusters,
             variant_hashes=args.variant_hashes,
             haploblock_hashes=args.haploblock_hashes,
             out=args.out)

    except Exception as e:
        # details on the issue should be in the exception name, print it to stderr and die
        sys.stderr.write("ERROR in " + script_name + " : " + repr(e) + "\n")
        sys.exit(1)