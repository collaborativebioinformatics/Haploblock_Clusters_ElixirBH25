import os
import logging
import pathlib


# set up logger, using inherited config, in case we get called as a module
logger = logging.getLogger(__name__)

import subprocess


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


def variant_counts_to_TSV(haploblock2count, out):
    """
    arguments:
    - haploblock2count: dict, key=(start, end), value=[mean, stdev]
    """
    with open(os.path.join(out, "variant_counts.tsv"), 'w') as f:
        # header
        f.write("START\tEND\tMEAN\tSTDEV\n")
        for (start, end) in haploblock2count:
            # mean and stdev in scientific notation with 3 significant digits
            mean = "{:.3g}".format(haploblock2count[(start, end)][0])
            stdev = "{:.3g}".format(haploblock2count[(start, end)][1])
            
            f.write(str(start) + "\t" + str(end) + "\t" + str(mean) + "\t" + str(stdev) + "\n")