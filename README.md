> Elixir BioHackathon November 3-7, 2025

# Haploblock_Clusters_ElixirBH25

This pipeline generates recombination-defined genomic hashes, unique identifiers for variant sets that include haplotype-resolved variant information. Every hash is a unique representation of an individual's genotype with a distinct combination of variants. We leverage large-scale population genomic data to train a model for a hash-based genotype-to-phenotype association prediction.


# How to use this repo

```
git clone https://github.com/collaborativebioinformatics/Haploblock_Clusters_ElixirBH25.git
cd Haploblock_Clusters_ElixirBH25/
```


# Data

All data listed below must be downloaded into `data/`:

```
cd data/
```

1. High-resolution recombination map from Halldorsson et al., 2019 with empirically defined genome-wide recombination rates:
```
## Download https://www.science.org/doi/suppl/10.1126/science.aau1043/suppl_file/aau1043_datas3.gz
## Upload the file to the data/ directory in this repo
gzip -d aau1043_datas3.gz
```

File `aau1043_datas3` contains averaged maternal and paternal recombination rates per genomic window.

2. 1000Genomes, HGSVC, Phase 3
For an example, we use chromosome 6. If you would like to run all chromosomes, download all of the phased VCF files.  
Phased VCF file of chromosome 6 (2548 samples):
```
wget https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000_genomes_project/release/20190312_biallelic_SNV_and_INDEL/ALL.chr6.shapeit2_integrated_snvindels_v2a_27022019.GRCh38.phased.vcf.gz
```

and index file (same point about chromosomes):
```
wget https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000_genomes_project/release/20190312_biallelic_SNV_and_INDEL/ALL.chr6.shapeit2_integrated_snvindels_v2a_27022019.GRCh38.phased.vcf.gz.tbi
```

3. Reference sequence for chromosome 6 (GRCh38), must be bgzipped:
```
wget https://hgdownload.soe.ucsc.edu/goldenPath/hg38/chromosomes/chr6.fa.gz
gzip -d chr6.fa.gz
## If you do not have bgzip, install it (see install_dependencies.txt)
bgzip chr6.fa
```

Optionally for testing, if you want to run the pipeline for one population from 1000Genomes, you will also need a TSV file with the list of samples in data/, eg. CHB with 113 samples downloaded from: https://www.internationalgenome.org/data-portal/population/CHB


# Run the pipeline

To create genomic hashes you can run the pipeline step-by-step or use a Docker container, as described below:


## Use Docker

Build the Docker image:
```
docker build -t haploblock-pipeline .
```

Run the Docker container interactively:
```
docker run -it --rm \
    -v host_mount_point/data:/app/Haploblock_Clusters_ElixirBH25/data \
    haploblock-pipeline
```

Once inside the container, the pipeline can be run with:
```
python3 main.py --config config/default.yaml [optional arguments]
```

The default configuration file example is provided in [config/default.yaml](haploblock_pipeline/config/default.yaml). This configuration file defines all input paths, chromosome settings, and parallelization parameters.

Optional pipeline arguments:
```
--threads <N> : Number of threads to use (default: the number of available CPUs minus 1)
--step <1|2|3|4|5> : Execute only a specific pipeline step (default: all steps)
```

Note: The pipeline can be also run inside a Docker container in the non-interactive mode with:
```
docker run --rm \
    -v host_mount_point/data:/app/Haploblock_Clusters_ElixirBH25/data \
    -v host_mount_point/config/default.yaml:/app/Haploblock_Clusters_ElixirBH25/config/default.yaml \
    haploblock-pipeline \
    python3 main.py --config config/default.yaml
```


## Run the pipeline step-by-step


### Configure Python environment and install dependencies

For the list of required dependencies and how to install them, please go to [requirements.md](requirements.md) and follow the instructions *carefully*.

The pipeline requires two Python libraries (numpy, pyyaml), as well as samtools, bcftools, htslib (see https://www.htslib.org/) and MMSeqs2 (https://github.com/soedinglab/MMseqs2). All must be simlinked in `/usr/bin` or exported to PATH.


#### 1. Generate haploblock boundaries and haploblock hashes for chromosome 6 using the Halldorsson2019 genetic map:

```
python3 haploblock_pipeline/step1_haploblocks.py \
    --recombination_file data/Halldorsson2019/aau1043_datas3 \
    --chr 6 \
    --out out_dir/
```

This step creates **haploblock boundaries**, a TSV file (with header) with 2 columns (START END), as well as **haploblock hashes**, a TSV file (with header) with 3 columns (START END HASH). Haploblock hashes are unique identifiers for haploblocks, i.e., strings of 1/0s of length equal to the number of haploblocks.

For a detailed description of how haplobock boundaries are defined see [https://github.com/jedrzejkubica/Elixir-BH-2025](https://github.com/jedrzejkubica/Elixir-BH-2025) as well as Halldorsson2019.


#### 2. Generate haploblock phased FASTA files and individual hashes:

This step uses the 1000Genomes phased VCF file to create a phased VCF with a subset of variants for every haploblock, then it generates a phased FASTA file for each individual.

Here is the instruction for **one haploblock (TNFa)** (example):
```
python haploblock_pipeline/step2_phased_sequences.py \
    --boundaries_file data/haploblock_boundaries_chr6_TNFa.tsv \
    --vcf data/ALL.chr6.shapeit2_integrated_snvindels_v2a_27022019.GRCh38.phased.vcf.gz \
    --ref data/chr6.fa.gz \
    --chr_map data/chr_map \
    --chr 6 \
    --variants data/variants_of_interest.txt \
    --out out_dir/TNFa/
```

Here is the instruction for **all haploblocks**:
```
python haploblock_pipeline/step2_phased_sequences.py \
    --boundaries_file out_dir/haploblock_boundaries_chr6.tsv \
    --vcf data/ALL.chr6.shapeit2_integrated_snvindels_v2a_27022019.GRCh38.phased.vcf.gz \
    --ref data/chr6.fa.gz \
    --chr_map data/chr_map \
    --chr 6 \
    --variants data/variants_of_interest.txt \
    --out out_dir/
```

NOTE: VCF file has "6" instead of "chr6", which is required by bcftools consensus, therefore create a file chr_map with one chromosome number-to-name mapping per line (e.g., "6 chr6") and provide it with --chr_map.

This script uses bcftools and bgzip to extract regions corresponding to haploblock boundaries (--boundaries_file) from a population VCF file (--vcf). Specify variants of interest in a file with one variant per line (--variants), they all must be in the same haploblock and in format: "chr(number only):position". We also calculate the mean and average of the number of variants per haploblock, they are saved in **out_dir/variant_counts.tsv** (with 4 columns: START, END, MEAN, STDEV). We assign individual hashes, ie integer numbers of lenght variants digits, each corresponding to variant of interest: 1 if variant in the sample or 0 otherwise, they are saved in **out_dir/variant_hashes.tsv** (with two columns: INDIVIDUAL HASH).

Optionally, if you want to run it for one population, use the TSV file with samples from 1000Genomes (data/igsr-chb.tsv.tsv) and run:
```
python haploblock_pipeline/step2_phased_sequences.py [...] --samples_file data/igsr-chb.tsv.tsv
```

This step creates a temporary directory in --out_dir, where it creates a consensus haploblock phased sequences for both haploids for each individual by applying common variants (min allele-frequence must be > 0.05) from previously generated phasesd VCF files to the reference sequence (--ref). 

Then generate one merged phased FASTA file per haploblock:

Here is the instruction for **one haploblock (TNFa)**:
```
haploblock_pipeline/step3_merge_fasta.py out_dir/TNFa/tmp/consensus_fasta out_dir/TNFa/haploblock_phased_seq_merged

## Remove the temporary directory out_dir/TNFa/tmp/ directory if you want
```

Here is the instruction for **all haploblocks**:
```
haploblock_pipeline/step3_merge_fasta.py out_dir/tmp/consensus_fasta out_dir/haploblock_phased_seq_merged

## Remove the temporary directory out_dir/tmp/ if you want
```

#### 3. Generate haploblock clusters:

This step generates haploblock clusters with MMSeqs2, make sure MMSeqs2 is available in PATH.

Here is the instruction for **one haploblock (TNFa)**:
```
python haploblock_pipeline/step4_clusters.py \
    --boundaries_file data/haploblock_boundaries_chr6_TNFa.tsv \
    --merged_consensus_dir out_dir/TNFa/haploblock_phased_seq_merged \
    --variant_counts out_dir/TNFa/variant_counts.tsv \
    --chr 6 \
    --out out_dir/TNFa/
```

Here is the instruction for **all haploblocks**:
```
python haploblock_pipeline/step4_clusters.py \
    --boundaries_file out_dir/haploblock_boundaries_chr6.tsv \
    --merged_consensus_dir out_dir/haploblock_phased_seq_merged \
    --variant_counts out_dir/variant_counts.tsv \
    --chr 6 \
    --out out_dir/
```

This step uses previously generated haploblock phased sequences (--merged_consensus_dir) and variant counts (--variant_counts), based on which it calculates MMSeqs2 parameters: min sequence identify and coverage fraction. For each haploblock it generates a cluster TSV file in directory **out_dir/clusters/**.


#### If you are not interested in SNPs, but only in haploblock clusters, you can stop here. Also, if you are interested in comparing these results to other groups of sequences, you can pull the cluster representatives from the merged fasta file (see step 2).

### 4. Generate variant hashes:

Here is the instruction for **one haploblock (TNFa)**:
```
python haploblock_pipeline/step4_clusters.py \
    --clusters out_dir/TNFa/clusters/chr6_31480875-31598421_cluster.tsv \
    --variant_hashes out_dir/TNFa/variant_hashes.tsv \
    --haploblock_hashes out_dir/haploblock_hashes_chr6.tsv \
    --chr 6 \
    --out out_dir/TNFa/
```

Here is the instruction for **all haploblocks**, please run it separately for each cluster TSV file:
```
python haploblock_pipeline/step4_clusters.py \
    --clusters out_dir/clusters/cluster_file.tsv \
    --variant_hashes out_dir/variant_hashes.tsv \
    --haploblock_hashes out_dir/haploblock_hashes_chr6.tsv \
    --chr 6 \
    --out out_dir/
```

This step generates variants hashes, 64-character strings of 0/1s. Each individual hash contains:
- strand hash: 4 characters
- chromosome hash: 10 characters
- haploblock hash: 20 characters
- cluster hash: 20 characters
- variant hash: the number of SNPs of interest

The output is a TSV file (with two columns: INDIVIDUAL HASH) in out_dir/individual_hashes.tsv.


## Example Model

Heres an example of how one could take individual hashes and feed them into a vector-based model with initial weighting...


# Results

## Testing the pipeline

We found 1399 haploblocks in chromosome 6. See [haploblock_boundaries_chr6.tsv](data/haploblock_boundaries_chr6.tsv) for these haploblock boundaries (high recombination rates defined as **rate > 10*average**).

We generated haploblock phased sequences (format: sample_chr_region_start-end_hap0/1.fa) for the CBH, PUR and GBR populations for the following genomic regions:
- 10 random haploblocks of chromosome 6
- 5 random haploblocks of chromosome 6
- the haploblocks overlapping with TNFa
- the haploblocks overlapping with genes related to human height

+ haploblock phased sequences and haploblock hashes for TNFa for all populations


# System requirements

During the hackathon we ran the pipeline on a Linux-based machine with 8 CPU cores and 16 GB of memory.


# Acknowlegdements

This work was supported by ELIXIR, the research infrastructure for life science data, and conducted at the ELIXIR BioHackathon Europe.


# References

1. Palsson, G., Hardarson, M.T., Jonsson, H. et al. Complete human recombination maps. Nature 639, 700–707 (2025). https://doi.org/10.1038/s41586-024-08450-5

2. Bjarni V. Halldorsson et al., Characterizing mutagenic effects of recombination through a sequence-level genetic map. Science363,eaau1043 (2019). DOI:10.1126/science.aau1043

3. Yu-Hsiang Tseng, Sumit Walia, Yatish Turakhia, "Ultrafast and ultralarge multiple sequence alignments using TWILIGHT", Bioinformatics, Volume 41, Issue Supplement_1, July 2025, Pages i332–i341, doi: 10.1093/bioinformatics/btaf212
