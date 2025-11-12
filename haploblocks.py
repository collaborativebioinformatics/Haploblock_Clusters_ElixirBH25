#!/usr/bin/env python3
import sys
import logging
import argparse
import pathlib
import numpy as np
from multiprocessing import Pool, cpu_count
import data_parser

logger = logging.getLogger(__name__)

HAPLOBLOCK_HASH_LENGTH = 20
#PARALLEL_THRESHOLD = 1_000_000  # Enable multiprocessing if > 1M haploblocks
PARALLEL_THRESHOLD = 1000  # DEV only option Trigger parallelization when >1000 haploblock

def _make_hash(i: int) -> str:
    """Helper function for multiprocessing."""
    return np.binary_repr(i, width=HAPLOBLOCK_HASH_LENGTH)


def generate_haploblock_hashes(haploblock_boundaries: list[tuple[int, int]]) -> dict[tuple[int, int], str]:
    """
    Generate unique binary hashes for each haploblock boundary.

    Uses multiprocessing for large datasets to speed up hash generation.

    Args:
        haploblock_boundaries: List of (start, end) tuples.

    Returns:
        Dict mapping (start, end) -> binary hash string.
    """
    n_blocks = len(haploblock_boundaries)
    logger.info(f"Generating hashes for {n_blocks:,} haploblocks")

    if n_blocks == 0:
        logger.warning("No haploblocks provided.")
        return {}

    # Decide execution strategy
    if n_blocks > PARALLEL_THRESHOLD:
        logger.info(f"Large dataset detected (> {PARALLEL_THRESHOLD:,} blocks). Using parallel generation.")
        with Pool(cpu_count()) as pool:
            hashes = pool.map(_make_hash, range(n_blocks))
    else:
        # Vectorized version for smaller datasets
        hashes = [np.binary_repr(i, width=HAPLOBLOCK_HASH_LENGTH) for i in range(n_blocks)]

    return dict(zip(haploblock_boundaries, hashes))


def haploblocks_to_tsv(haploblock_boundaries: list[tuple[int, int]], chrom: str, out_dir: pathlib.Path) -> None:
    """Save haploblock boundaries to TSV file."""
    output_file = out_dir / f"haploblock_boundaries_chr{chrom}.tsv"
    with output_file.open('w') as f:
        f.write("START\tEND\n")
        f.writelines(f"{start}\t{end}\n" for start, end in haploblock_boundaries)


def haploblock_hashes_to_tsv(haploblock2hash: dict[tuple[int, int], str], chrom: str, out_dir: pathlib.Path) -> None:
    """Save haploblock hashes to TSV file."""
    output_file = out_dir / f"haploblock_hashes_chr{chrom}.tsv"
    with output_file.open('w') as f:
        f.write("START\tEND\tHASH\n")
        f.writelines(f"{start}\t{end}\t{hash_}\n" for (start, end), hash_ in haploblock2hash.items())


def main(recombination_file: pathlib.Path, chrom: str, out_dir: pathlib.Path) -> None:
    logger.info(f"Parsing recombination file: {recombination_file}")
    haploblock_boundaries = data_parser.parse_recombination_rates(recombination_file, chrom)

    # Ensure output directory exists
    out_dir.mkdir(parents=True, exist_ok=True)

    logger.info("Generating haploblock hashes")
    haploblock2hash = generate_haploblock_hashes(haploblock_boundaries)

    logger.info("Saving haploblock boundaries and hashes")
    haploblocks_to_tsv(haploblock_boundaries, chrom, out_dir)
    haploblock_hashes_to_tsv(haploblock2hash, chrom, out_dir)


# -------------------------------------------------------------------------
# CUDA Perspective (Optional, Commented Out)
# -------------------------------------------------------------------------
# For extremely large datasets and GPU-enabled environments, this
# section shows how hashing could be performed on CUDA for massive speedups.
#
# from numba import cuda
#
# @cuda.jit
# def generate_hashes_cuda(output, width):
#     i = cuda.grid(1)
#     if i < output.size:
#         val = i
#         for j in range(width):
#             bit = (val >> (width - j - 1)) & 1
#             output[i, j] = bit
#
# def generate_haploblock_hashes_cuda(n_blocks):
#     import numpy as np
#     output = np.zeros((n_blocks, HAPLOBLOCK_HASH_LENGTH), dtype=np.uint8)
#     threads_per_block = 256
#     blocks_per_grid = (n_blocks + (threads_per_block - 1)) // threads_per_block
#     generate_hashes_cuda[blocks_per_grid, threads_per_block](output, HAPLOBLOCK_HASH_LENGTH)
#     cuda.synchronize()
#     return ["".join(str(bit) for bit in row) for row in output]
# -------------------------------------------------------------------------


if __name__ == "__main__":
    script_name = pathlib.Path(sys.argv[0]).name

    logging.basicConfig(
        format="%(asctime)s %(levelname)s %(name)s: %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
        level=logging.DEBUG,
    )

    parser = argparse.ArgumentParser(
        prog=script_name,
        description="Generate haploblock boundaries and hashes from recombination data."
    )
    parser.add_argument('--recombination_file', type=pathlib.Path, required=True,
                        help='Path to recombination file (Halldorsson et al., 2019)')
    parser.add_argument('--chr', type=str, required=True, help='Chromosome name or number')
    parser.add_argument('--out', type=pathlib.Path, required=True, help='Output folder path')

    args = parser.parse_args()

    try:
        main(args.recombination_file, args.chr, args.out)
    except Exception as e:
        sys.stderr.write(f"ERROR in {script_name}: {repr(e)}\n")
        sys.exit(1)

