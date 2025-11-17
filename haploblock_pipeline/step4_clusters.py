#!/usr/bin/env python3
import os
import sys
import logging
import pathlib
import subprocess
from concurrent.futures import ThreadPoolExecutor, as_completed

import data_parser

logger = logging.getLogger(__name__)

# ----------------------------------------------------------------------
# MMseqs2 params
# ----------------------------------------------------------------------
def calculate_mmseq_params(variant_counts_file: pathlib.Path):
    haploblock2min_id = {}
    haploblock2cov_fraction = {}

    with open(variant_counts_file, "r") as f:
        header = f.readline()
        if not header.startswith("START\t"):
            raise ValueError(f"Variant counts file missing header: {header.strip()}")

        for line in f:
            start, end, mean, stdev = line.strip().split("\t")
            hap_len = int(end) - int(start)
            haploblock2min_id[(start, end)] = 1 - (float(mean) / hap_len)
            haploblock2cov_fraction[(start, end)] = 1 - (682 / hap_len)

    return haploblock2min_id, haploblock2cov_fraction

# ----------------------------------------------------------------------
# Run clustering per FASTA
# ----------------------------------------------------------------------

def compute_clusters(input_fasta: str, out: str, min_seq_id: float, cov_fraction: float, cov_mode: int,
                     chrom: str, start: str, end: str):

    # Resolve EVERYTHING to absolute paths
    input_fasta = str(pathlib.Path(input_fasta).resolve())

    output_prefix = pathlib.Path(out) / "clusters" / f"chr{chrom}_{start}-{end}"
    output_prefix = str(output_prefix.resolve())

    tmp_dir = pathlib.Path(out) / "tmp"
    tmp_dir = str(tmp_dir.resolve())

    cmd = [
        "mmseqs", "easy-cluster",
        input_fasta,
        output_prefix,
        tmp_dir,
        "--min-seq-id", str(min_seq_id),
        "-c", str(cov_fraction),
        "--cov-mode", str(cov_mode),
        "--remove-tmp-files", "1"
    ]

    logger.debug("Running: %s", " ".join(cmd))
    subprocess.run(cmd, check=True)


# ----------------------------------------------------------------------
# Main workflow
# ----------------------------------------------------------------------
def run_clusters(boundaries_file: pathlib.Path, merged_consensus_dir: pathlib.Path,
                 variant_counts_file: pathlib.Path, chrom: str, out_dir: pathlib.Path,
                 cov_mode: int, threads: int = None):

    out_dir.mkdir(parents=True, exist_ok=True)
    (out_dir / "clusters").mkdir(exist_ok=True)
    (out_dir / "tmp").mkdir(exist_ok=True)

    haploblock_boundaries = data_parser.parse_haploblock_boundaries(boundaries_file)
    haploblock2min_id, haploblock2cov_fraction = calculate_mmseq_params(variant_counts_file)

    logger.info("Found %d haploblocks", len(haploblock_boundaries))
    logger.info("Computing clusters with MMseqs2 using %s threads...", threads or "auto")

    max_workers = threads or max(1, os.cpu_count() - 1)
    futures = []

    with ThreadPoolExecutor(max_workers=max_workers) as executor:
        for start, end in haploblock_boundaries:
            input_fasta = merged_consensus_dir / f"chr{chrom}_region_{start}-{end}.fa"
            futures.append(
                executor.submit(
                    compute_clusters,
                    str(input_fasta),
                    str(out_dir),
                    haploblock2min_id[(start, end)],
                    haploblock2cov_fraction[(start, end)],
                    cov_mode,
                    chrom, start, end
                )
            )

        for fut in as_completed(futures):
            fut.result()  # propagate exceptions if any

    logger.info("All clusters computed successfully.")

# ----------------------------------------------------------------------
# Pipeline wrapper
# ----------------------------------------------------------------------
def run(boundaries_file, merged_consensus_dir, variant_counts, chr, out, cov_mode=0, threads=None):
    run_clusters(
        pathlib.Path(boundaries_file),
        pathlib.Path(merged_consensus_dir),
        pathlib.Path(variant_counts),
        str(chr),
        pathlib.Path(out),
        cov_mode=cov_mode,
        threads=threads
    )

# ----------------------------------------------------------------------
# CLI
# ----------------------------------------------------------------------
if __name__ == "__main__":
    logging.basicConfig(
        format="%(asctime)s %(levelname)s %(name)s: %(message)s",
        level=logging.INFO
    )
    parser = argparse.ArgumentParser(description="Cluster haploblock consensus sequences using MMseqs2")
    parser.add_argument("--boundaries_file", type=pathlib.Path, required=True)
    parser.add_argument("--merged_consensus_dir", type=pathlib.Path, required=True)
    parser.add_argument("--variant_counts", type=pathlib.Path, required=True)
    parser.add_argument("--chr", type=str, required=True)
    parser.add_argument("--out", type=pathlib.Path, required=True)
    parser.add_argument("--cov_mode", type=int, default=0)
    parser.add_argument("--threads", type=int, default=None)

    args = parser.parse_args()

    run(
        boundaries_file=args.boundaries_file,
        merged_consensus_dir=args.merged_consensus_dir,
        variant_counts=args.variant_counts,
        chr=args.chr,
        out=args.out,
        cov_mode=args.cov_mode,
        threads=args.threads
    )
