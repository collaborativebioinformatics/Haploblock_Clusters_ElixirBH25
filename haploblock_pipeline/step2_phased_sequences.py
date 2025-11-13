#!/usr/bin/env python3
import sys
import logging
import argparse
import pathlib
import subprocess
import numpy as np
import data_parser
from multiprocessing import Pool, cpu_count
from typing import Tuple, List, Dict

logger = logging.getLogger(__name__)


def count_variants(vcf: pathlib.Path) -> Tuple[int, int]:
    """
    Count the number of variants in a bgzipped VCF file.

    Returns:
        (count_hap0, count_hap1)
    """
    result = subprocess.run(
        ["bcftools", "query", "-f", "[ %GT]", str(vcf)],
        capture_output=True,
        text=True,
        check=True
    )

    count_0 = count_1 = 0
    for token in result.stdout.split():
        if token in {".", "./.", ".|."}:
            continue

        if "|" in token:
            left, right = token.split("|")
        elif "/" in token:
            left, right = token.split("/")
        else:
            continue

        count_0 += left == "1"
        count_1 += right == "1"

    return count_0, count_1


def generate_consensus_fasta(ref_fasta: pathlib.Path,
                             vcf: pathlib.Path,
                             out_dir: pathlib.Path) -> Tuple[pathlib.Path, pathlib.Path]:
    """
    Apply variants from VCF to the reference to create phased consensus FASTA files.
    Uses bcftools consensus (external command) — keep it sequential per call, but we
    can run multiple calls in parallel from a Pool.
    """
    tmp_dir = out_dir / "tmp" / "consensus_fasta"
    tmp_dir.mkdir(parents=True, exist_ok=True)

    # Ensure the vcf name doesn't carry a trailing ".vcf" portion from repeated stems
    vcf_name = vcf.stem
    if vcf_name.endswith(".vcf"):
        vcf_name = vcf_name[:-4]

    output_hap0 = tmp_dir / f"{vcf_name}_hap0.fa"
    output_hap1 = tmp_dir / f"{vcf_name}_hap1.fa"

    for hap, outfile in [(1, output_hap0), (2, output_hap1)]:
        with outfile.open("w") as f_out:
            subprocess.run(
                ["bcftools", "consensus", "-H", str(hap), "-f", str(ref_fasta), str(vcf)],
                stdout=f_out,
                check=True
            )

    return output_hap0, output_hap1


def generate_variant_hashes(variants: List[str],
                            vcf: pathlib.Path,
                            chrom: str,
                            haploblock_boundaries,
                            samples: List[str]) -> Dict[str, str]:
    """
    Generate binary variant presence hashes for all samples and haplotypes.

    This function queries the VCF once (per region) for all samples and
    builds the variant2hash mapping. It's left sequential since bcftools
    handles many-sample queries efficiently. If you need to parallelize
    across many regions, call this function from multiple processes.
    """
    first_variant_pos = variants[0]
    start = end = None
    for (bound_start, bound_end) in haploblock_boundaries:
        if int(bound_start) <= int(first_variant_pos) <= int(bound_end):
            start, end = bound_start, bound_end
            break

    if start is None:
        raise ValueError(f"Variant {first_variant_pos} not found in any haploblock")

    variant_indices = {str(pos): i for i, pos in enumerate(variants)}
    variant2hash = {
        f"{sample}_chr{chrom}_region_{start}-{end}_hap{h}": ["0"] * len(variants)
        for sample in samples for h in (0, 1)
    }

    query_region = f"{chrom}:{first_variant_pos}-{end}"
    result = subprocess.run(
        ["bcftools", "query", "-f", "%CHROM\t%POS[\t%GT]\n",
         "-s", ",".join(samples), "--force-samples", "-r", query_region, str(vcf)],
        capture_output=True, text=True, check=True
    )

    for line in result.stdout.splitlines():
        chrom_col, pos, *genotypes = line.split("\t")
        if pos not in variant_indices:
            continue
        idx = variant_indices[pos]

        for sample_idx, genotype in enumerate(genotypes):
            if "|" not in genotype:
                continue
            hap0, hap1 = genotype.split("|")
            sample = samples[sample_idx]
            if hap0 == "1":
                variant2hash[f"{sample}_chr{chrom}_region_{start}-{end}_hap0"][idx] = "1"
            if hap1 == "1":
                variant2hash[f"{sample}_chr{chrom}_region_{start}-{end}_hap1"][idx] = "1"

    # convert hash lists to strings
    return {k: "".join(v) for k, v in variant2hash.items()}

    # -------------------------------------------------------------------------
    # Perspective: CUDA Acceleration
    # -------------------------------------------------------------------------
    # For very large datasets (many haploblocks × samples), the variant-hash
    # generation step can be offloaded to a GPU. This allows thousands of 
    # haplotype/variant comparisons to occur simultaneously.
    #
    # Example concept (requires CuPy or PyTorch):
    #
    # import cupy as cp
    #
    # def generate_variant_hashes_cuda(variants, genotypes):
    #     # Convert genotypes (numpy array of 0/1 values) to GPU array
    #     geno_gpu = cp.array(genotypes, dtype=cp.uint8)
    #     # Compute hashes in parallel
    #     hash_gpu = cp.packbits(geno_gpu, axis=1)
    #     # Transfer back to CPU
    #     return cp.asnumpy(hash_gpu)
    #
    # In practice:
    # - Each sample × haplotype combination can be represented as a binary vector.
    # - The GPU computes all binary encodings and bitwise operations at once.
    #
    # Similarly, for consensus FASTA generation:
    # - The reference sequence and VCF variants can be loaded as GPU tensors.
    # - Base substitution and mask operations (A→T, G→C, etc.) are trivial for CUDA.
    #
    # CUDA acceleration could reduce runtime from hours → minutes on large cohorts.
    # -------------------------------------------------------------------------


def write_tsv_variant_hashes(variant2hash: Dict[str, str], out_dir: pathlib.Path) -> None:
    """Write variant hashes to TSV."""
    out_path = out_dir / "variant_hashes.tsv"
    with out_path.open("w") as f:
        f.write("INDIVIDUAL\tHASH\n")
        f.writelines(f"{ind}\t{hash_}\n" for ind, hash_ in variant2hash.items())


def write_tsv_variant_counts(haploblock2count: Dict[tuple, Tuple[float, float]], out_dir: pathlib.Path) -> None:
    """Write mean and stdev variant counts to TSV."""
    out_path = out_dir / "variant_counts.tsv"
    with out_path.open("w") as f:
        f.write("START\tEND\tMEAN\tSTDEV\n")
        for (start, end), (mean, stdev) in haploblock2count.items():
            f.write(f"{start}\t{end}\t{mean:.3g}\t{stdev:.3g}\n")


# --- Worker for per-sample work (to be run in Pool) -----------------------
def _process_sample_for_region(args):
    """
    Worker function to:
      - extract sample VCF for a given region
      - count variants for both haplotypes
      - generate consensus FASTAs for the sample
    Returns: (sample, count0, count1)
    """
    (sample, region_vcf, region_fasta, out_dir) = args
    try:
        sample_vcf = data_parser.extract_sample_from_vcf(region_vcf, sample, out_dir)
        count_0, count_1 = count_variants(sample_vcf)
        # Generate consensus FASTA (writes to out_dir/tmp/consensus_fasta)
        generate_consensus_fasta(region_fasta, sample_vcf, out_dir)
        return (sample, count_0, count_1)
    except Exception as e:
        logger.exception("Error processing sample %s: %s", sample, e)
        # Return zeros on failure to avoid crashing the whole pool (caller can detect)
        return (sample, 0, 0)


def run_phased_sequences(boundaries_file: pathlib.Path,
         samples_file: pathlib.Path | None,
         vcf: pathlib.Path,
         ref: pathlib.Path,
         chr_map: pathlib.Path,
         chrom: str,
         variants_file: pathlib.Path,
         out_dir: pathlib.Path,
         workers: int = None) -> None:

    # Ensure output dirs exist
    tmp_dir = out_dir / "tmp"
    if tmp_dir.exists():
        logger.error(f"Output directory {out_dir} already exists; remove it first.")
        raise FileExistsError(f"{out_dir} exists")
    (tmp_dir / "consensus_fasta").mkdir(parents=True, exist_ok=True)

    logger.info("Parsing haploblock boundaries")
    haploblock_boundaries = data_parser.parse_haploblock_boundaries(boundaries_file)
    logger.info("Found %d haploblocks", len(haploblock_boundaries))

    logger.info("Parsing samples")
    samples = (data_parser.parse_samples(samples_file)
               if samples_file else data_parser.parse_samples_from_vcf(vcf))
    logger.info("Found %d samples", len(samples))

    logger.info("Parsing variants of interest")
    variants = data_parser.parse_variants_of_interest(variants_file)

    # Generate the variant2hash for the region containing the first variant (kept sequential)
    variant2hash = generate_variant_hashes(variants, vcf, chrom, haploblock_boundaries, samples)
    write_tsv_variant_hashes(variant2hash, out_dir)

    # Decide on workers
    cpu_cores = cpu_count()
    if workers is None or workers <= 0:
        workers = cpu_cores
    logger.info("Using %d worker(s) for per-sample processing (available cores: %d)", workers, cpu_cores)

    # Per-haploblock processing: extract region vcf/fasta, then process samples (parallelizable)
    haploblock2count = {}
    for (start, end) in haploblock_boundaries:
        logger.info("Processing haploblock %s-%s", start, end)
        region_vcf = data_parser.extract_region_from_vcf(vcf, chrom, chr_map, start, end, out_dir)
        region_fasta = data_parser.extract_region_from_fasta(ref, chrom, start, end, out_dir)

        # Prepare work items for samples
        work_items = [(sample, region_vcf, region_fasta, out_dir) for sample in samples]

        if workers == 1:
            # Sequential processing (fallback)
            counts = []
            for item in work_items:
                sample, c0, c1 = _process_sample_for_region(item)
                counts.extend([c0, c1])
        else:
            # Parallel processing across samples (Pool)
            with Pool(workers) as pool:
                results = pool.map(_process_sample_for_region, work_items)
            counts = []
            for (_, c0, c1) in results:
                counts.extend([c0, c1])

        # Record mean and stdev for this haploblock
        haploblock2count[(start, end)] = [np.mean(counts) if counts else 0.0,
                                          np.std(counts) if counts else 0.0]

    write_tsv_variant_counts(haploblock2count, out_dir)


if __name__ == "__main__":
    logging.basicConfig(
        format="%(asctime)s %(levelname)s %(name)s: %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
        level=logging.DEBUG,
    )
    logger = logging.getLogger(pathlib.Path(sys.argv[0]).name)

    parser = argparse.ArgumentParser(description="Generate haploblock phased sequences and hashes.")
    parser.add_argument("--boundaries_file", type=pathlib.Path, required=True)
    parser.add_argument("--vcf", type=pathlib.Path, required=True)
    parser.add_argument("--ref", type=pathlib.Path, required=True)
    parser.add_argument("--chr_map", type=pathlib.Path, required=True)
    parser.add_argument("--chr", type=str, required=True)
    parser.add_argument("--variants", type=pathlib.Path, required=True)
    parser.add_argument("--out", type=pathlib.Path, required=True)
    parser.add_argument("--samples_file", type=pathlib.Path, required=False)
    parser.add_argument("--workers", type=int, required=False, default=None,
                        help="Number of worker processes for per-sample parallelism (default: cpu_count()).")

    args = parser.parse_args()
    try:
        run_phased_sequences(args.boundaries_file, args.samples_file, args.vcf, args.ref,
             args.chr_map, args.chr, args.variants, args.out, args.workers)
    except Exception as e:
        sys.stderr.write(f"ERROR in {pathlib.Path(sys.argv[0]).name}: {repr(e)}\n")
        sys.exit(1)

