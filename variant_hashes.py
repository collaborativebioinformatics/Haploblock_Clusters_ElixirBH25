import os
import sys
import logging
import argparse
import pathlib


# set up logger, using inherited config, in case we get called as a module
logger = logging.getLogger(__name__)

def parse_clusters(clusters_file):
    """
    Parses clusters file from MMSeqs2 (no header) with 2 columns: representative, member
    assing unique ids for each cluster.
    We want to create unique cluster IDs based on cluster representatives,
    and match members to cluster IDs.

    arguments:
    - clusters_file

    returns:
    - member2cluster: dict, key=member, value=unique clusterID
    - num_clusters: int
    """
    try:
        f = open(clusters_file, 'r')
    except Exception as e:
        logger.error("Opening provided clusters file %s: %s", clusters_file, e)
        raise Exception("Cannot open provided clusters file")
    
    representative2cluster = {}
    member2representative = {}
    member2cluster = {}
    num_clusters = 0
    for line in f:
        split_line = line.rstrip().split('\t')

        if len(split_line) != 2:
            logger.error("Clusters file %s has bad line (not 2 tab-separated fields): %s",
                         clusters_file, line)
            raise Exception("Bad line in the clusters file")
        
        (representative, member) = split_line
        
        if not representative in representative2cluster:
            representative2cluster[representative] = num_clusters
            num_clusters += 1

        member2representative[member] = representative
    
    for member in member2representative:
        representative = member2representative[member]
        cluster = representative2cluster[representative]
        member2cluster[member] = cluster

    return(member2cluster, num_clusters)


def parse_individual_hashes():
    pass

def parse_haploblock_hashes():
    pass

def generate_variant_hases():
    # generate cluster hashes
    pass



def main(clusters_file, individual_hashes, haploblock_hashes):

    logger.info("Parsing clusters")
    (member2cluster, num_clusters) = parse_clusters(clusters_file)
    logger.info("Found %i clusters with %i members in total", num_clusters, len(member2cluster))


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
    parser.add_argument('--individual_hashes',
                        help='Path to file with individual hashes',
                        type=pathlib.Path,
                        required=True)
    parser.add_argument('--haploblock_hashes',
                        help='Path to file with haploblock hashes',
                        type=pathlib.Path,
                        required=True)

    args = parser.parse_args()

    try:
        main(clusters_file=args.clusters,
             individual_hashes=args.individual_hashes,
             haploblock_hashes=args.haploblock_hashes)

    except Exception as e:
        # details on the issue should be in the exception name, print it to stderr and die
        sys.stderr.write("ERROR in " + script_name + " : " + repr(e) + "\n")
        sys.exit(1)