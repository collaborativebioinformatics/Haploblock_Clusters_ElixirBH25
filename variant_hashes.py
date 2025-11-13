#!/usr/bin/env python3
import os
import sys
import logging
import argparse
import pathlib
import numpy as np
from concurrent.futures import ThreadPoolExecutor, as_completed

import data_parser

logger = logging.getLogger(__name__)

CLUSTER_HASH_LENGTH = 20


def generate_cluster_hashes(clusters):
    """
    Generate a unique binary hash for each cluster.
    Parallelization is not needed here — small overhead.
    """
    cluster2hash = {cluster: np.binary_repr(idx, width=CLUSTER_HASH_LENGTH)
                    for idx, cluster in enumerate(clusters)}
    return cluster2hash


def generate_individual_hash(individual, individual2cluster, cluster2hash,
                             variant2hash, haploblock2hash, chr_hash):
    """Generate hash string for a single individual."""
    strand = individual[-1]
    strand_hash = "0001" if strand == "0" else "0010"

    # Extract start-end robustly
    individual_split = individual.split("_")
    region_str = individual_split[3].replace(".fa", "").replace(".fasta", "").replace(".vcf", "")
    start, end = region_str.split("-")
    haploblock_hash = haploblock2hash[(start, end)]
    cluster_hash = cluster2hash[individual2cluster[individual]]
    variant_hash = variant2hash[individual]

    return individual, strand_hash + chr_hash + haploblock_hash + cluster_hash + variant_hash


def generate_individual_hashes_parallel(individual2cluster, cluster2hash, variant2hash,
                                        haploblock2hash, chr_hash, max_workers=None):
    """
    Parallelized version: generate hashes for all individuals using ThreadPoolExecutor.
    """
    haploblock2hash = {(str(s), str(e)): h for (s, e), h in haploblock2hash.items()}
    individual2hash = {}

    max_workers = max_workers or (os.cpu_count() - 1 or 1)
    futures = []

    with ThreadPoolExecutor(max_workers=max_workers) as executor:
        for individual in individual2cluster:
            futures.append(
                executor.submit(
                    generate_individual_hash,
                    individual, individual2cluster, cluster2hash,
                    variant2hash, haploblock2hash, chr_hash
                )
            )
        for fut in as_completed(futures):
            ind, h = fut.result()
            individual2hash[ind] = h

    return individual2hash


def hashes_to_TSV(individual2hash, out):
    output_path = os.path.join(out, "individual_hashes.tsv")
    with open(output_path, 'w') as f:
        f.write("INDIVIDUAL\tHASH\n")
        for individual, hash_val in individual2hash.items():
            f.write(f"{individual}\t{hash_val}\n")
    logger.info(f"Individual hashes written to {output_path}")


def main(clusters_file, variant_hashes, haploblock_hashes, chr, out):
    logger.info("Parsing clusters")
    individual2cluster, clusters = data_parser.parse_clusters(clusters_file)
    logger.info(f"Found {len(clusters)} clusters with {len(individual2cluster)} individuals")

    logger.info("Parsing variant hashes")
    variant2hash = data_parser.parse_variant_hashes(variant_hashes)

    logger.info("Parsing haploblock hashes")
    haploblock2hash = data_parser.parse_haploblock_hashes(haploblock_hashes)

    logger.info("Generating individual hashes (parallelized)")
    chr_hash = np.binary_repr(int(chr))
    cluster2hash = generate_cluster_hashes(clusters)
    individual2hash = generate_individual_hashes_parallel(
        individual2cluster, cluster2hash, variant2hash, haploblock2hash, chr_hash
    )

    logger.info("Saving individual hashes")
    hashes_to_TSV(individual2hash, out)


# ----------------------------------------------------------------------
# CUDA perspective
# ----------------------------------------------------------------------
#    Bitwise hash computation is perfect for GPU acceleration:
#    - Represent hashes as uint8/uint32 arrays on GPU
#    - Use CuPy or PyTorch for batch concatenation of strand, chr, haploblock, cluster, variant hashes
#    - Can scale to millions of individuals efficiently
#    - Final conversion to string can be done on CPU before TSV export


if __name__ == "__main__":
    script_name = os.path.basename(sys.argv[0])
    logging.basicConfig(
        format="%(asctime)s %(levelname)s %(name)s: %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
        level=logging.INFO,
    )
    logger = logging.getLogger(script_name)

    parser = argparse.ArgumentParser(
        prog=script_name, description="Generate individual variant hashes per haploblock"
    )
    parser.add_argument('--clusters', required=True, type=pathlib.Path,
                        help='Path to clusters file generated with MMSeqs2')
    parser.add_argument('--variant_hashes', required=True, type=pathlib.Path,
                        help='Path to file with variant hashes')
    parser.add_argument('--haploblock_hashes', required=True, type=pathlib.Path,
                        help='Path to file with haploblock hashes')
    parser.add_argument('--chr', required=True, type=str, help='Chromosome number')
    parser.add_argument('--out', required=True, type=pathlib.Path,
                        help='Output folder path')

    args = parser.parse_args()

    try:
        main(
            clusters_file=args.clusters,
            variant_hashes=args.variant_hashes,
            haploblock_hashes=args.haploblock_hashes,
            chr=args.chr,
            out=args.out,
        )
    except Exception as e:
        sys.stderr.write(f"ERROR in {script_name}: {repr(e)}\n")
        sys.exit(1)

