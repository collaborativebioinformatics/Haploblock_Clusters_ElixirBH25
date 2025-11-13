#!/usr/bin/env python3
import argparse
import logging
import sys
import yaml
from pathlib import Path

from haploblocks_pipeline import (
    step1_haploblocks,
    step2_phased_sequences,
    step3_merge_fasta,
    step4_clusters,
    step5_variant_hashes,
)
from haploblocks_pipeline.utils.logging import setup_logger


def load_config(config_path):
    """Load YAML configuration file."""
    with open(config_path, "r") as f:
        return yaml.safe_load(f)


def main():
    parser = argparse.ArgumentParser(description="Run the Haploblocks Pipeline")
    parser.add_argument("--config", type=Path, help="Path to YAML configuration file")

    # Optional CLI overrides
    parser.add_argument("--step", type=str, help="Pipeline step (1–5 or all)")
    parser.add_argument("--threads", type=int, help="Number of CPU threads")

    args = parser.parse_args()
    logger = setup_logger()

    # Load config
    if args.config:
        logger.info(f"Loading configuration from {args.config}")
        cfg = load_config(args.config)
    else:
        logger.error("No config file provided. Use --config config/default.yaml")
        sys.exit(1)

    # Allow CLI overrides
    step = args.step or cfg["pipeline"]["step"]
    threads = args.threads if args.threads else cfg["pipeline"]["threads"]
    if threads == "auto":
        import os
        threads = max(1, (os.cpu_count() or 2) - 1)

    logger.info(f"Starting Haploblocks pipeline (Step: {step}, Threads: {threads})")

    try:
        # Step 1
        if step in ["1", "all"]:
            step1_haploblocks.run(
                recombination_file=cfg["data"]["recombination_file"],
                chr=cfg["chromosome"]["number"],
                out=cfg["outputs"]["out_dir"],
                threads=threads
            )

        # Step 2
        if step in ["2", "all"]:
            step2_phased_sequences.run(
                boundaries_file=cfg["data"]["boundaries_file"],
                vcf=cfg["data"]["vcf"],
                ref=cfg["data"]["ref"],
                chr_map=cfg["data"]["chr_map"],
                chr=cfg["chromosome"]["number"],
                variants=cfg["data"]["variants"],
                out=Path(cfg["outputs"]["out_dir"]) / "TNFa",
                threads=threads
            )

        # Step 3
        if step in ["3", "all"]:
            step3_merge_fasta.run(
                input_dir=Path(cfg["outputs"]["out_dir"]) / "TNFa" / "tmp" / "consensus_fasta",
                output_dir=Path(cfg["outputs"]["out_dir"]) / "TNFa" / "haploblock_phased_seq_merged",
                threads=threads
            )

        # Step 4
        if step in ["4", "all"]:
            step4_clusters.run(
                boundaries_file=cfg["data"]["boundaries_file"],
                merged_consensus_dir=cfg["outputs"]["merged_consensus_dir"],
                variant_counts=cfg["outputs"]["variant_counts"],
                chr=cfg["chromosome"]["number"],
                out=Path(cfg["outputs"]["out_dir"]) / "TNFa",
                threads=threads
            )

        # Step 5
        if step in ["5", "all"]:
            step5_variant_hashes.run(
                clusters=cfg["outputs"]["clusters"],
                variant_hashes=cfg["outputs"]["variant_hashes"],
                haploblock_hashes=cfg["outputs"]["haploblock_hashes"],
                chr=cfg["chromosome"]["number"],
                out=Path(cfg["outputs"]["out_dir"]) / "TNFa",
                threads=threads
            )

    except Exception as e:
        logger.error(f"Pipeline failed: {e}")
        sys.exit(1)

    logger.info("Pipeline finished successfully!")


if __name__ == "__main__":
    main()

