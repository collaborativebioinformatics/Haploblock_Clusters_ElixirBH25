import os
import sys
import logging
import argparse
import pathlib
import subprocess

import numpy

import data_parser

# set up logger, using inherited config, in case we get called as a module
logger = logging.getLogger(__name__)


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


def generate_consensus_fasta(fasta, vcf, out):
    """
    Apply variants from VCF to reference sequence

    Generates the following files in out/ :
    - ref_chr6.fa.gz.fai
    - ref_chr6.fa.gz.gzi
    - chr{chr}_region_{start}-{end}.fa.gz
    """
    output_hap0 = os.path.join(out, "tmp", pathlib.Path(vcf.stem).stem + "_hap0.fa")  # removes .vcf.gz
    output_hap1 = os.path.join(out, "tmp", pathlib.Path(vcf.stem).stem + "_hap1.fa")  # removes .vcf.gz
    
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


def generate_variant_hashes(variants, vcf, chr, haploblock_boundaries, samples):
    """
    Generate a "variant hash" for every variant of interest,
    ie an integer number with len(variants) digits, each corresponding to variant of interes:
    1 if variant in vcf, 0 otherwise

    the way we do it here is probably not the most efficient

    arguments:
    - vcf: pathlib.Path to bgzipped vcf
    - variants: list of variants of interest

    returns:
    - individual2variantHash: dict, key=sampleRegionHap, value=int(variant hash)
    """        
    # find haploblock with variants of interest
    first_variant = variants[0]  # we assume that the user provided variants within the same haploblock
    start = 0
    end = 0
    for (bound_start, bound_end) in haploblock_boundaries:
        if int(bound_start) <= int(first_variant) <= int(bound_end):
            start = bound_start
            end = bound_end
            break
        else:
            logger.error("Cannot find variant %s", first_variant)
            raise Exception("Variants not in any haploblock")
    
    # initialize all variant hashes to string of "0"s
    individual2variantHash = {}
    for sample in samples:
        individual_1 = f"{sample}_chr{chr}_region_{start}-{end}_hap1"
        individual_2 = f"{sample}_chr{chr}_region_{start}-{end}_hap2"

        individual2variantHash[individual_1] = "0" * len(variants)
        individual2variantHash[individual_2] = "0" * len(variants)

    # search for genotypes in VCF
    query = chr + ":" + ",".join(variants)
    sub_vcf = subprocess.run(["bcftools", "query",
                             "-f", "%CHROM\t%POS[\t%GT]\n",
                             "-s", ",".join(samples),
                             "--force-samples",
                             "-r", query,
                             vcf],
                             capture_output=True,
                             text=True,
                             check=True).stdout.splitlines()
    
    # assign variant hashes
    for (variant_count, line) in enumerate(sub_vcf):
        line_split = line.split("\t")
        chrom, pos, *genotypes = line_split
        for (sample_count, genotype) in enumerate(genotypes):
            sample = samples[sample_count]
            if '|' not in genotype:
                logger.warning("variant %s:%s is unphased or missing, skipping it",  chrom, pos)
                continue
            hap1, hap2 = genotype.split('|')

            individual_1 = f"{sample}_chr{chr}_region_{start}-{end}_hap1"
            individual_2 = f"{sample}_chr{chr}_region_{start}-{end}_hap2"

            if hap1 == "1":
                print("hap1 = 1")
                individual2variantHash[individual_1][variant_count] = "1"
            if hap2 == "1":
                print("hap2 = 2")
                individual2variantHash[individual_2][variant_count] = "1"
    
    return(individual2variantHash)


def variant_hashes_to_TSV(individual2variantHash, out):
    """
    arguments:
    - : dict, key=sampleRegionHap, value=variant hash
    """
    with open(os.path.join(out, "variant_hashes.tsv"), 'w') as f:
        # header
        f.write("INDIVIDUAL\tVAR_HASH\n")
        for individual in individual2variantHash:
            f.write(individual + "\t" + individual2variantHash[individual] + "\n")


def main(boundaries_file, samples_file, vcf, ref, chr_map, chr, out):
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
    
    # create out/ and a tmp/ directory in out/ for temporary files
    if os.path.exists(os.path.join(out, "tmp")):
        logger.error(f"Output directory {os.path.join(out)} exists, please remove it")
        raise Exception("Output directory exists")
    os.makedirs(out)
    if os.path.exists(os.path.join(out, "tmp")):
        logger.info(f"Temporary directory {os.path.join(out, 'tmp')} exists, removing it")
        subprocess.run(["rm", "-r", os.path.join(out, "tmp")],
                       check=True)
    os.mkdir(os.path.join(out, "tmp"))

    logger.info("Parsing haploblock boundaries")
    haploblock_boundaries = data_parser.parse_haploblock_boundaries(boundaries_file)
    logger.info("Found %i haploblocks", len(haploblock_boundaries))

    logger.info("Parsing samples")
    samples = data_parser.parse_samples(samples_file)
    # samples = data_parser.parse_samples_from_vcf(vcf)
    logger.info("Found %i samples", len(samples))

    # dict for variant counts, key=(start, end), value=list(mean, stdev)
    haploblock2count = {}

    # we assume that the user provided variants within the same haploblock,
    # need to check it later
    variants = ["31480875"]
    individual2variantHash = generate_variant_hashes(variants, vcf, chr, haploblock_boundaries, samples)
    variant_hashes_to_TSV(individual2variantHash, out)
    exit(1)

    for (start, end) in haploblock_boundaries:
        logger.info(f"Generating phased VCF for haploblock {start}-{end}")
        region_vcf = data_parser.extract_region_from_vcf(vcf, chr, chr_map, start, end, out)

        logger.info(f"Generating phased fasta for haploblock {start}-{end}")
        region_fasta = data_parser.extract_region_from_fasta(ref, chr, start, end, out)

        # list of the number of variants per sample and per haplotype
        haploblock_counts = []

        logger.info(f"Generating consensus fasta files for haploblock {start}-{end}")
        for sample in samples:
            # logger.info(f"Generating phased VCF for haploblock {start}-{end} for sample %s", sample)
            sample_vcf = data_parser.extract_sample_from_vcf(region_vcf, sample, out)

            # calculate the number variants in the sample
            (count_0, count_1) = count_variants(sample_vcf)
            haploblock_counts.append(count_0)
            haploblock_counts.append(count_1)

            (sample_hap1, sample_hap2) = generate_consensus_fasta(region_fasta, sample_vcf, out)

        # calculate mean and stdev of the number variants
        mean = sum(haploblock_counts) / len(haploblock_counts)
        stdev = numpy.std(haploblock_counts)

        haploblock2count[(start, end)] = [mean, stdev]
    
    data_parser.variant_counts_to_TSV(haploblock2count, out)


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

    args = parser.parse_args()

    try:
        main(boundaries_file=args.boundaries_file,
             samples_file=args.samples_file,
             vcf=args.vcf,
             ref=args.ref,
             chr_map=args.chr_map,
             chr=args.chr,
             out=args.out)

    except Exception as e:
        # details on the issue should be in the exception name, print it to stderr and die
        sys.stderr.write("ERROR in " + script_name + " : " + repr(e) + "\n")
        sys.exit(1)