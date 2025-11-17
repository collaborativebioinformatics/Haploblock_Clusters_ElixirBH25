import os
import logging
import pathlib
import subprocess
import numpy as np

# set up logger, using inherited config
logger = logging.getLogger(__name__)


# -------------------------------------------------------------------------
# Recombination rates (vectorized)
# -------------------------------------------------------------------------
def parse_recombination_rates(recombination_file, chromosome):
    """
    Parses recombination rates from Halldorsson et al., 2019
    and returns haploblock boundaries as list of tuples (start, end).
    """
    if not str(chromosome).startswith("chr"):
        chromosome = f"chr{chromosome}"

    try:
        # load file, skip comments
        data = np.genfromtxt(
            recombination_file,
            delimiter="\t",
            dtype=None,
            names=True,
            encoding=None,
            comments="#",
        )
    except Exception as e:
        logger.error("Cannot open recombination file %s: %s", recombination_file, e)
        raise

    # filter chromosome
    chr_mask = data['Chr'] == chromosome
    chr_data = data[chr_mask]

    if len(chr_data) == 0:
        logger.warning("No positions found for %s in %s", chromosome, recombination_file)
        return []

    positions = chr_data['Begin'].astype(int)
    rates = chr_data['cMperMb'].astype(float)

    avg_rate = np.mean(rates)
    high_rate_mask = rates > 10 * avg_rate
    high_rate_positions = positions[high_rate_mask]

    if len(high_rate_positions) == 0:
        return [(1, positions[-1])]

    haploblocks = [(1, high_rate_positions[0])]
    haploblocks += [(high_rate_positions[i - 1], high_rate_positions[i]) for i in range(1, len(high_rate_positions))]
    haploblocks.append((high_rate_positions[-1], positions[-1]))

    logger.info("Found %d haploblocks (NumPy optimized)", len(haploblocks))
    return haploblocks


# -------------------------------------------------------------------------
# Haploblock boundaries
# -------------------------------------------------------------------------
def parse_haploblock_boundaries(boundaries_file):
    haploblocks = []
    with open(boundaries_file, 'r') as f:
        header = f.readline()
        if not header.startswith("START\t"):
            raise Exception(f"Boundaries file {boundaries_file} missing header")
        for line in f:
            start, end = line.rstrip().split('\t')
            haploblocks.append((start, end))
    return haploblocks


# -------------------------------------------------------------------------
# Samples parsing
# -------------------------------------------------------------------------
def parse_samples(samples_file):
    samples = []
    with open(samples_file, 'r') as f:
        header = f.readline()
        if not header.startswith("Sample name\t"):
            raise Exception(f"Samples file {samples_file} missing header")
        for line in f:
            sample, *rest = line.rstrip().split('\t')
            if not (sample.startswith("HG") or sample.startswith("NA")):
                raise Exception(f"Invalid sample line: {line}")
            samples.append(sample)
    return samples


def parse_samples_from_vcf(vcf):
    result = subprocess.run(
        ["bcftools", "query", "-l", vcf],
        check=True,
        capture_output=True,
        text=True,
    )
    return result.stdout.splitlines()


# -------------------------------------------------------------------------
# Variants of interest
# -------------------------------------------------------------------------
def parse_variants_of_interest(variants_file):
    variants = []
    with open(variants_file, 'r') as f:
        for line in f:
            parts = line.rstrip().split(":")
            if len(parts) != 2:
                raise Exception(f"Bad line in variants file: {line}")
            _, pos = parts
            variants.append(pos)
    return variants


# -------------------------------------------------------------------------
# Extract region/sample
# -------------------------------------------------------------------------
def extract_region_from_vcf(vcf, chr, chr_map, start, end, out):
    if chr.startswith("chr"):
        chr = chr.replace("chr", "")

    tmp_dir = pathlib.Path(out) / "tmp"
    tmp_dir.mkdir(parents=True, exist_ok=True)

    temporary_vcf = tmp_dir / f"{chr}_region_{start}-{end}.vcf.gz"

    subprocess.run([
        "bcftools", "view",
        "-r", f"{chr}:{start}-{end}",
        "--min-af", "0.05",
        str(vcf),
        "-o", str(temporary_vcf)
    ], check=True)

    subprocess.run(["bcftools", "index", str(temporary_vcf)], check=True)

    output_vcf = tmp_dir / f"chr{chr}_region_{start}-{end}.vcf"
    output_index = tmp_dir / f"chr{chr}_region_{start}-{end}.vcf.gz.csi"

    subprocess.run([
        "bcftools", "annotate",
        "--rename-chrs", str(chr_map),
        str(temporary_vcf)
    ], stdout=open(output_vcf, "w"), check=True)

    subprocess.run(["bgzip", str(output_vcf)], check=True)
    output_vcf_bgzip = pathlib.Path(str(output_vcf) + ".gz")
    subprocess.run(["bcftools", "index", "-c", "-o", str(output_index), str(output_vcf_bgzip)], check=True)

    return output_vcf_bgzip


def extract_sample_from_vcf(vcf, sample, out):
    tmp_dir = pathlib.Path(out) / "tmp"
    tmp_dir.mkdir(parents=True, exist_ok=True)
    output_vcf = tmp_dir / f"{sample}_{vcf.stem}.gz"

    subprocess.run([
        "bcftools", "view",
        "--force-samples",
        "-s", sample,
        "-o", str(output_vcf),
        str(vcf)
    ], check=True)

    subprocess.run(["bcftools", "index", str(output_vcf)], check=True)
    return output_vcf


def extract_region_from_fasta(fasta, chr, start, end, out):
    tmp_dir = pathlib.Path(out) / "tmp"
    tmp_dir.mkdir(parents=True, exist_ok=True)
    subprocess.run(["samtools", "faidx", str(fasta)], check=True)
    output_fasta = tmp_dir / f"chr{chr}_region_{start}-{end}.fa"
    subprocess.run(["samtools", "faidx", str(fasta), f"chr{chr}:{start}-{end}"],
                   stdout=open(output_fasta, "w"), check=True)
    return output_fasta


# -------------------------------------------------------------------------
# Clusters
# -------------------------------------------------------------------------
def parse_clusters(clusters_file):
    individual2cluster = {}
    clusters = []
    rep2cluster = {}
    num_clusters = 0

    with open(clusters_file, 'r') as f:
        for line in f:
            rep, indiv = line.rstrip().split('\t')
            if rep not in rep2cluster:
                rep2cluster[rep] = num_clusters
                clusters.append(num_clusters)
                num_clusters += 1
            individual2cluster[indiv] = rep2cluster[rep]

    return individual2cluster, clusters


# -------------------------------------------------------------------------
# Variant / Haploblock hashes
# -------------------------------------------------------------------------
def parse_variant_hashes(variant_hashes):
    d = {}
    with open(variant_hashes, 'r') as f:
        header = f.readline()
        for line in f:
            indiv, h = line.rstrip().split('\t')
            d[indiv] = h
    return d


def parse_haploblock_hashes(haploblock_hashes):
    d = {}
    with open(haploblock_hashes, 'r') as f:
        header = f.readline()
        for line in f:
            start, end, h = line.rstrip().split('\t')
            d[(start, end)] = h
    return d

