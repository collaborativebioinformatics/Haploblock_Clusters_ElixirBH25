import os
import sys
import logging
import argparse
import pathlib

import numpy


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

    Generates the following files in out/ :
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
    temporary_vcf = os.path.join(out, f"{chr}_region_{start}-{end}.vcf.gz")

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
    output_vcf = os.path.join(out, f"chr{chr}_region_{start}-{end}.vcf")
    output_index = os.path.join(out, f"chr{chr}_region_{start}-{end}.vcf.gz.csi")

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

    output_vcf = os.path.join(out, sample + "_" + vcf.stem + ".gz")

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


def count_variants(vcf):
    """
    Count the number of variants in the vcf file
    """
    # extract genotype (GT) strings
    GTs = subprocess.run(["bcftools", "query",
                          "-f", "[ %GT]",
                          vcf],
                          capture_output=True,
                          text=True,
                          check=True)

    count_0 = 0  # left haplotype
    count_1 = 0  # right haplotype

    for line in GTs.stdout.splitlines():
        line = line.strip()
        if not line or line == "." or line == "./." or line == ".|.":
            continue  # skip missing genotypes

        # Handle phased and unphased cases
        if "|" in line:
            parts = line.split("|")
        elif "/" in line:
            parts = line.split("/")
        else:
            continue  # skip malformed entries

        if len(parts) != 2:
            continue  # unexpected format, skip

        (left, right) = parts

        if left == "1":
            count_0 += 1
        if right == "1":
            count_1 += 1

    return(count_0, count_1)


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
    
    output_fasta = os.path.join(out, f"chr{chr}_region_{start}-{end}.fa")
    # extract region start-end from reference fasta
    subprocess.run(["samtools", "faidx",
                    fasta,
                    f"chr{chr}:{start}-{end}"],
                    stdout=open(output_fasta, "w"),
                    check=True)

    return(output_fasta)


def generate_consensus_fasta(fasta, vcf, out):
    """
    Apply variants from VCF to reference sequence

    Generates the following files in out/ :
    - ref_chr6.fa.gz.fai
    - ref_chr6.fa.gz.gzi
    - chr{chr}_region_{start}-{end}.fa.gz
    """
    output_hap1 = os.path.join(out, pathlib.Path(vcf.stem).stem + "_hap1.fa")  # removes .vcf.gz
    output_hap2 = os.path.join(out, pathlib.Path(vcf.stem).stem + "_hap2.fa")  # removes .vcf.gz
    
    # create a consensus sequence (fasta) from reference and variants extracted from VCF
    # haploid sequence 1
    subprocess.run(["bcftools", "consensus",
                    "-H", "1",
                    "-f", fasta,
                    vcf],
                    stdout=open(output_hap1, "w"),
                    check=True)
    
    # haploid sequence 2
    subprocess.run(["bcftools", "consensus",
                    "-H", "2",
                    "-f", fasta,
                    vcf],
                    stdout=open(output_hap2, "w"),
                    check=True)

    return(output_hap1, output_hap2)


def save_haploblock_variant_counts(haploblock2count, out_file):
    """
    arguments:
    - haploblock2count: dict, key=(start, end), value=[mean, stdev]
    """
    with open(out_file, 'w') as f:
        # header
        f.write("START\tEND\tMEAN\tSTDEV\n")
        for (start, end) in haploblock2count:
            # mean and stdev in scientific notation with 3 significant digits
            mean = "{:.3g}".format(haploblock2count[(start, end)][0])
            stdev = "{:.3g}".format(haploblock2count[(start, end)][1])
            
            f.write(str(start) + "\t" + str(end) + "\t" + str(mean) + "\t" + str(stdev) + "\n")


def main(boundaries_file, samples_file, vcf, ref, chr_map, chr, out, variant_counts_file):
    # sanity check
    if not os.path.exists(boundaries_file):
        logger.error(f"File {boundaries_file} does not exist.")
        raise Exception("File does not exist")
    if not os.path.exists(samples_file):
        logger.error(f"File {samples_file} does not exist.")
        raise Exception("File does not exist")
    if not os.path.exists(vcf):
        logger.error(f"File {vcf} does not exist.")
        raise Exception("File does not exist")
    if not os.path.exists(ref):
        logger.error(f"File {ref} does not exist.")
        raise Exception("File does not exist")
    if not os.path.exists(chr_map):
        logger.error(f"File {chr_map} does not exist.")
        raise Exception("File does not exist")

    logger.info("Parsing haploblock boundaries")
    haploblock_boundaries = parse_haploblock_boundaries(boundaries_file)
    logger.info("Found %i haploblock boundaries", len(haploblock_boundaries))

    logger.info("Parsing samples")
    samples = parse_samples(samples_file)
    # samples = parse_samples_from_vcf(vcf)
    logger.info("Found %i samples", len(samples))

    # dict for variant counts, key=(start, end), value=list(mean, stdev)
    haploblock2count = {}

    for (start, end) in haploblock_boundaries:
        logger.info(f"Generating phased VCF for haploblock {start}-{end}")
        region_vcf = extract_region_from_vcf(vcf, chr, chr_map, start, end, out)

        logger.info(f"Generating phased fasta for haploblock {start}-{end}")
        region_fasta = extract_region_from_fasta(ref, chr, start, end, out)

        # list for the number variants (separate for each haplotype)
        haploblock_counts = []

        logger.info(f"Generating consensus fasta files for haploblock {start}-{end}")
        for sample in samples:
            # logger.info(f"Generating phased VCF for haploblock {start}-{end} for sample %s", sample)
            sample_vcf = extract_sample_from_vcf(region_vcf, sample, out)

            # calculate of the number variants
            (count_0, count_1) = count_variants(sample_vcf)
            haploblock_counts.append(count_0)
            haploblock_counts.append(count_1)

            (sample_hap1, sample_hap2) = generate_consensus_fasta(region_fasta, sample_vcf, out)

        # calculate mean and stdev of the number variants
        mean = sum(haploblock_counts) / len(haploblock_counts)
        stdev = numpy.std(haploblock_counts)

        haploblock2count[(start, end)] = [mean, stdev]
    
    save_haploblock_variant_counts(haploblock2count, variant_counts_file)


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
    
    parser.add_argument('--boundaries_file',
                        help='Path to boundaries file generated from Halldorsson et al., 2019',
                        type=pathlib.Path,
                        required=True)
    parser.add_argument('--samples_file',
                        help='Path to samples file from 1000Genomes',
                        type=pathlib.Path,
                        required=True)
    parser.add_argument('--vcf',
                        help='Path to phased VCF file (bgzipped) from 1000Genomes',
                        type=pathlib.Path,
                        required=True)
    parser.add_argument('--ref',
                        help='Path to reference sequence (bgzipped)',
                        type=pathlib.Path,
                        required=True)
    parser.add_argument('--chr_map',
                        help='Path to chr_map: one mapping per line, ie "6 chr6"',
                        type=pathlib.Path,
                        required=True)
    parser.add_argument('--chr',
                        help='chromosome',
                        type=str,
                        required=True)
    parser.add_argument('--out',
                        help='Path to output folder',
                        type=pathlib.Path,
                        required=True)
    parser.add_argument('--variant_counts_file',
                        help='Path to a file to variant counts (mean, stdev)',
                        type=pathlib.Path,
                        required=True)

    args = parser.parse_args()

    try:
        main(boundaries_file=args.boundaries_file,
             samples_file=args.samples_file,
             vcf=args.vcf,
             ref=args.ref,
             chr_map=args.chr_map,
             chr=args.chr,
             out=args.out,
             variant_counts_file=args.variant_counts_file)

    except Exception as e:
        # details on the issue should be in the exception name, print it to stderr and die
        sys.stderr.write("ERROR in " + script_name + " : " + repr(e) + "\n")
        sys.exit(1)