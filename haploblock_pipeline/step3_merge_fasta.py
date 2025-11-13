#!/usr/bin/env python3
import os
import sys
import logging
import pathlib
import subprocess
from haploblocks_pipeline.utils.logging import setup_logger

logger = setup_logger()

def run(input_dir: pathlib.Path, output_dir: pathlib.Path, threads: int | None = None, clean: bool = False):
    """
    Wrapper to call the merge_fasta_per_region.sh script with proper threads handling.
    """
    if threads is None or threads <= 0:
        threads = max(1, (os.cpu_count() or 2) - 1)

    if not input_dir.exists():
        logger.error("Input directory does not exist: %s", input_dir)
        raise FileNotFoundError(input_dir)

    output_dir.mkdir(parents=True, exist_ok=True)

    cmd = [
        "bash",
        str(pathlib.Path(__file__).parent / "merge_fasta_per_region.sh"),
        str(input_dir),
        str(output_dir),
    ]
    if clean:
        cmd.append("--clean")
    # Append threads as last argument
    cmd.append(str(threads))

    logger.info("Running merge FASTA per region with %d threads...", threads)
    try:
        subprocess.run(cmd, check=True)
    except subprocess.CalledProcessError as e:
        logger.error("FASTA merge failed: %s", e)
        sys.exit(1)

    logger.info("FASTA merge completed successfully.")

