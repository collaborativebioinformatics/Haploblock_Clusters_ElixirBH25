import os
import logging
import pathlib

import subprocess

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
    
    # add first haploblock
    haploblock_boundaries.append((1, high_rates_positions[0]))
    
    for i in range(1, len(high_rates)):
        start = high_rates_positions[i - 1]
        end = high_rates_positions[i]
        haploblock_boundaries.append((start, end))
    
    # add last haploblock
    haploblock_boundaries.append((high_rates_positions[-1], positions[-1]))

    logger.info("Found %i haploblocks", len(haploblock_boundaries))
    
    return(haploblock_boundaries)


def parse_haploblock_boundaries(boundaries_file):
    """
    Parses haploblock boundaries with 2 columns: start end

    arguments:
    - boundaries_file
    
    returns:
    - haploblock_boundaries: list of tuples with haploblock boundaries (start, end)
    """
    haploblock_boundaries = []

    try:
        f = open(boundaries_file, 'r')
    except Exception as e:
        logger.error("Opening provided boundaries file %s: %s", boundaries_file, e)
        raise Exception("Cannot open provided boundaries file")
    
    # skip header
    line = f.readline()
    if not line.startswith("START\t"):
        logging.error("boundaries file %s is headerless? expecting headers but got %s",
                      boundaries_file, line)
        raise Exception("boundaries file problem")
    
    for line in f:
        split_line = line.rstrip().split('\t')

        if len(split_line) != 2:
            logger.error("Boundaries file %s has bad line (not 2 tab-separated fields): %s",
                         boundaries_file, line)
            raise Exception("Bad line in the boundaries file")
        
        (start, end) = split_line

        haploblock_boundaries.append((start, end))

    return(haploblock_boundaries)


def parse_samples(samples_file):
    """
    Parse samples file for a population from 1000Genomes

    arguments:
    - samples file with 9 columns:
    Sample name, Sex, Biosample ID, Population code, Population name, Superpopulation code,
    Superpopulation name, Population elastic ID, Data collections

    returns:
    - list of sample names
    """
    samples = []

    try:
        f = open(samples_file, 'r')
    except Exception as e:
        logger.error("Opening provided samples file %s: %s", samples_file, e)
        raise Exception("Cannot open provided samples file")
    
    # skip header
    line = f.readline()
    if not line.startswith("Sample name\t"):
        logging.error("samples file %s is headerless? expecting headers but got %s",
                      samples_file, line)
        raise Exception("samples file problem")
    
    for line in f:
        split_line = line.rstrip().split('\t')
        
        (sample, *_) = split_line

        if not (sample.startswith("HG") or sample.startswith("NA")):
            logger.error("Samples file %s has bad line (not 9 tab-separated fields): %s",
                         samples_file, line)
            raise Exception("Bad line in the samples file")

        samples.append(sample)

    return(samples)


def parse_samples_from_vcf(vcf):
    """
    Parse samples file for a population VCF from 1000Genomes

    arguments:
    - vcf file

    returns:
    - list of sample names
    """
    samples = []

    try:
        f = open(vcf, 'r')
    except Exception as e:
        logger.error("Opening provided vcf file %s: %s", vcf, e)
        raise Exception("Cannot open provided vcf file")
    
    samples = subprocess.run(["bcftools", "query",
                              "-l",
                              vcf],
                              check=True,
                              capture_output=True,
                              text=True).stdout.splitlines()

    return(samples)


def parse_variants_of_interest(variants_file):
    """
    Parses variants of interest file with one variant per line in format chr:pos

    arguments:
    - variants_file
    
    returns:
    - variants: list of string variant positions
    """
    variants = []

    try:
        f = open(variants_file, 'r')
    except Exception as e:
        logger.error("Opening provided variants file %s: %s", variants_file, e)
        raise Exception("Cannot open provided variants file")
    
    for line in f:
        line_split = line.rstrip().split(":")

        if len(line_split) != 2:
            logger.error("Variants file %s has bad line: %s",
                         variants_file, line)
            raise Exception("Bad line in the variants file")
        
        chr, position = line_split

        variants.append(position)

    return(variants)


def extract_region_from_vcf(vcf, chr, chr_map, start, end, out):
    """
    Extract a specific region from a VCF file

    Generates the following files in out/tmp/:
    - {chr}_region_{start}-{end}.vcf.gz
    - {chr}_region_{start}-{end}.vcf.gz.csi
    - chr{chr}_region_{start}-{end}.vcf
    - chr{chr}_region_{start}-{end}.vcf.gz
    - chr{chr}_region_{start}-{end}.vcf.gz.csi

    returns:
    - output_vcf: pathlib.Path to bgzipped vcf
    """
    if chr.startswith("chr"):
        chr = chr.replace("chr", "")

    # extract region start-end from VCF and index
    temporary_vcf = os.path.join(out, "tmp", f"{chr}_region_{start}-{end}.vcf.gz")

    subprocess.run(["bcftools", "view",
                    "-r", f"{chr}:{start}-{end}",
                    "--min-af", "0.05",
                    vcf,
                    "-o", temporary_vcf],
                    check=True)

    subprocess.run(["bcftools", "index",
                    temporary_vcf],
                    check=True)
    
    # VCF has 6 instead of chr6, which is required by bcftools consensus
    # create file chr_map: "6 chr6" one mapping per line
    # map chr6 to 6, bgzip and index
    output_vcf = os.path.join(out, "tmp", f"chr{chr}_region_{start}-{end}.vcf")
    output_index = os.path.join(out, "tmp", f"chr{chr}_region_{start}-{end}.vcf.gz.csi")

    subprocess.run(["bcftools", "annotate",
                    "--rename-chrs", chr_map,
                    temporary_vcf],
                    stdout=open(output_vcf, "w"),
                    check=True)
    
    subprocess.run(["bgzip",
                    output_vcf],
                    check=True)

    output_vcf_bgzip = output_vcf + ".gz"

    subprocess.run(["bcftools", "index",
                    "-c",
                    "-o", output_index,
                    output_vcf_bgzip],
                    check=True)
    
    return(pathlib.Path(output_vcf_bgzip))


def extract_sample_from_vcf(vcf, sample, out):
    """
    Extract a specific sample from a VCF file
    
    returns:
    - output_vcf: pathlib.Path to bgzipped VCF
    """

    output_vcf = os.path.join(out, "tmp", sample + "_" + vcf.stem + ".gz")

    # extract sample from VCF and index
    subprocess.run(["bcftools", "view",
                    "--force-samples",  # only warn about unknown subset samples
                    "-s", sample,
                    "-o", output_vcf,
                    vcf],
                    check=True)

    subprocess.run(["bcftools", "index",
                    output_vcf],
                    check=True)
    
    return(pathlib.Path(output_vcf))


def extract_region_from_fasta(fasta, chr, start, end, out):
    """
    Extract a specific region from a fasta file
    
    returns:
    - output_fasta: pathlib.Path to fasta with region start-end
    """
    # index reference
    subprocess.run(["samtools", "faidx",
                    fasta],
                    check=True)
    
    output_fasta = os.path.join(out, "tmp", f"chr{chr}_region_{start}-{end}.fa")
    # extract region start-end from reference fasta
    subprocess.run(["samtools", "faidx",
                    fasta,
                    f"chr{chr}:{start}-{end}"],
                    stdout=open(output_fasta, "w"),
                    check=True)

    return(output_fasta)


def parse_clusters(clusters_file):
    """
    Parses clusters file from MMSeqs2 (no header) with 2 columns: representative, individual
    assing unique ids for each cluster.
    We want to create unique cluster ID (starting at 0) based on cluster representatives,
    and match individual to cluster IDs.

    arguments:
    - clusters_file

    returns:
    - individual2cluster: dict, key=individual, value=unique clusterID
    - clusters: list of cluster IDs
    """
    try:
        f = open(clusters_file, 'r')
    except Exception as e:
        logger.error("Opening provided clusters file %s: %s", clusters_file, e)
        raise Exception("Cannot open provided clusters file")
    
    representative2cluster = {}
    individual2representative = {}
    individual2cluster = {}
    clusters = []
    num_clusters = 0
    for line in f:
        split_line = line.rstrip().split('\t')

        if len(split_line) != 2:
            logger.error("Clusters file %s has bad line (not 2 tab-separated fields): %s",
                         clusters_file, line)
            raise Exception("Bad line in the clusters file")
        
        (representative, individual) = split_line
        
        if not representative in representative2cluster:
            representative2cluster[representative] = num_clusters
            clusters.append(num_clusters)
            num_clusters += 1

        individual2representative[individual] = representative
    
    for individual in individual2representative:
        representative = individual2representative[individual]
        cluster = representative2cluster[representative]
        individual2cluster[individual] = cluster

    return(individual2cluster, clusters)


def parse_variant_hashes(variant_hashes):
    """
    Parses file with variant hashes with two columns: INDIVIDUAL HASH

    arguments:
    - variant_hashes: pathlib.Path
    
    returns:
    - variant2hash: dict, key=individual, key=hash
    """
    variant2hash = {}

    try:
        f = open(variant_hashes, 'r')
    except Exception as e:
        logger.error("Opening provided variant hashes file %s: %s", variant_hashes, e)
        raise Exception("Cannot open provided variant hashes file")
    
    # skip header
    line = f.readline()
    if not line.startswith("INDIVIDUAL\t"):
        logging.error("hashes file %s is headerless? expecting headers but got %s",
                      variant_hashes, line)
        raise Exception("hashes file problem")
    
    for line in f:
        split_line = line.rstrip().split('\t')

        if len(split_line) != 2:
            logger.error("Hashes file %s has bad line (not 2 tab-separated fields): %s",
                         variant_hashes, line)
            raise Exception("Bad line in the hashes file")
        
        (individual, hash) = split_line

        variant2hash[individual] = hash

    return(variant2hash)


def parse_haploblock_hashes(haploblock_hashes):
    """
    Parses file with haploblock hashes with 3 columns: START END HASH

    arguments:
    - haploblock_hashes: pathlib.Path
    
    returns:
    - haploblock2hash: dict, key=(start, end), key=hash
    """
    haploblock2hash = {}

    try:
        f = open(haploblock_hashes, 'r')
    except Exception as e:
        logger.error("Opening provided haploblock hashes file %s: %s", haploblock_hashes, e)
        raise Exception("Cannot open provided haploblock hashes file")
    
    # skip header
    line = f.readline()
    if not line.startswith("START\t"):
        logging.error("hashes file %s is headerless? expecting headers but got %s",
                      haploblock_hashes, line)
        raise Exception("hashes file problem")
    
    for line in f:
        split_line = line.rstrip().split('\t')

        if len(split_line) != 3:
            logger.error("Hashes file %s has bad line (not 3 tab-separated fields): %s",
                         haploblock_hashes, line)
            raise Exception("Bad line in the hashes file")
        
        (start, end, hash) = split_line

        haploblock2hash[(start, end)] = hash

    return(haploblock2hash)
