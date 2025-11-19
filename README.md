> Elixir BioHackathon November 3-7, 2025

# Haploblock_Clusters_ElixirBH25

# How to use this repo

```
git clone https://github.com/collaborativebioinformatics/Haploblock_Clusters_ElixirBH25.git
cd Haploblock_Clusters_ElixirBH25/
```


### Configure Python environment

Install via [Python venv](https://docs.python.org/3/library/venv.html) with the following command:

```
python3 -m venv --system-site-packages ~/pyEnv_ElixirBH2025
source ~/pyEnv_ElixirBH2025/bin/activate
pip install --upgrade pip
pip install numpy

```


### Install other dependencies

Please go to [install_dependencies.txt](install_dependencies.txt) and follow the instructions *carefully*.

Install samtools, bcftools, htslib (https://www.htslib.org/) and MMSeqs2 (https://github.com/soedinglab/MMseqs2). All must be simlinked in `/usr/bin` or exported to PATH.


# Data

All data listed below must be downloaded into `data/`:

```
cd data/
```

1. high-resolution recombination map from Halldorsson et al., 2019 with empirically defined recombination rates:
```
## Download https://www.science.org/doi/suppl/10.1126/science.aau1043/suppl_file/aau1043_datas3.gz
## Upload the file to the data/ directory in this repo
gzip -d aau1043_datas3.gz
```

File `aau1043_datas3` contains averaged maternal and paternal recombination rates.

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

Optionally, this is for testing, if you want to run this for one population from 1000Genomes, you will also need a TSV file with the list of samples in data/, eg. CHB (113 samples): https://www.internationalgenome.org/data-portal/population/CHB


# Workflow

TBD


# Run pipeline

```
## If you are in data/
cd ..
```


#### 1. Generate haploblock boundaries and hashes for chr6 using the Halldorsson2019 recombination map:

```
python haploblocks.py \
    --recombination_file data/Halldorsson2019/aau1043_datas3 \
    --chr 6 \
    --out out_dir/
```

This script creates **haploblock boundaries**, a TSV file (with header) with 2 columns (START END), as well as **haploblock hashes**, a TSV file (with header) with 3 columns (START END HASH). Haploblock hashes are unique identifiers for haploblocks, ie integers with len(haploblock_boundaries) digits 1/0s. All files are saved into `out_dir/`.

For a detailed description of how haplobock boundaries are defined see: [https://github.com/jedrzejkubica/Elixir-BH-2025](https://github.com/jedrzejkubica/Elixir-BH-2025) and Halldorsson2019.


#### 2. Generate haploblock phased fasta files (1000Genomes phased VCF -> Haploblock phased VCFs -> Phased fasta files) and individual hashes:

Here is the instruction for **one haploblock (TNFa)**:
```
python haploblock_phased_sequences.py \
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
python haploblock_phased_sequences.py \
    --boundaries_file out_dir/haploblock_boundaries_chr6.tsv \
    --vcf data/ALL.chr6.shapeit2_integrated_snvindels_v2a_27022019.GRCh38.phased.vcf.gz \
    --ref data/chr6.fa.gz \
    --chr_map data/chr_map \
    --chr 6 \
    --variants data/variants_of_interest.txt \
    --out out_dir/
```

NOTE: VCF file has "6" instead of "chr6", which is required by bcftools consensus, create file chr_map with one mapping per line (e.g., "6 chr6") and provide it using --chr_map.

Optionally, if you want to run it for one population, use the TSV file with samples from 1000Genomes (data/igsr-chb.tsv.tsv) and provide it as follows:
```
python haploblock_phased_sequences.py [...] --samples_file data/igsr-chb.tsv.tsv
```

This script uses bcftools and bgzip to extract regions corresponding to haploblock boundaries (--boundaries_file) from a population VCF file (--vcf). Specify variants of interest in a file with one variant per line (--variants), they all must be in the same haploblock and in format: "chr(number only):position". We also calculate the mean and average of the number of variants per haploblock, they are saved in **out_dir/variant_counts.tsv** (with 4 columns: START, END, MEAN, STDEV). We assign individual hashes, ie integer numbers of lenght variants digits, each corresponding to variant of interest: 1 if variant in the sample or 0 otherwise, they are saved in **out_dir/variant_hashes.tsv** (with two columns: INDIVIDUAL HASH)

It creates out_dir/tmp/ and generates consensus haploblock phased sequences for both haploids of each sample (e.g., `NA18531_chr6_region_711055-761032_hap1.fa`) by applying common variants (bcftools view `--min-af 0.05`) from previously generated VCF to reference sequence (--ref). 

Generate one merged phased fasta file per haploblock:

Here is the instruction for **one haploblock (TNFa)**:
```
./merge_fasta_per_region.sh out_dir/TNFa/tmp/consensus_fasta out_dir/TNFa/haploblock_phased_seq_merged

## Remove the out_dir/TNFa/tmp/ directory if you want
```

Here is the instruction for **all haploblocks**:
```
./merge_fasta_per_region.sh out_dir/tmp/consensus_fasta out_dir/haploblock_phased_seq_merged

## Remove the out_dir/tmp/ directory if you want
```

#### 3. Generate haploblock clusters with MMSeqs2:

Make sure MMSeqs is available in PATH.

Here is the instruction for **one haploblock (TNFa)**:
```
python clusters.py \
    --boundaries_file data/haploblock_boundaries_chr6_TNFa.tsv \
    --merged_consensus_dir out_dir/TNFa/haploblock_phased_seq_merged \
    --variant_counts out_dir/TNFa/variant_counts.tsv \
    --chr 6 \
    --out out_dir/TNFa/
```

Here is the instruction for **all haploblocks**:
```
python clusters.py \
    --boundaries_file out_dir/haploblock_boundaries_chr6.tsv \
    --merged_consensus_dir out_dir/haploblock_phased_seq_merged \
    --variant_counts out_dir/variant_counts.tsv \
    --chr 6 \
    --out out_dir/
```

This uses previously generated haploblock phased sequences (--merged_consensus_dir) and variant counts (--variant_counts), based on which it calculates MMSeqs parameters: min sequence identify and coverage fraction. For each haploblock it generates a cluster TSV file in directory **out_dir/clusters/**.

### 3a. If you are not interested in SNPs, and only background haploblock clusters, you can stop here.

#### 3a1. If you are interested in comparing these results to other groups of sequences, you can pull the cluster representatives from the merged fasta file (see step 2).

### 4. Generate variant hashes

Here is the instruction for **one haploblock (TNFa)**:
```
python variant_hashes.py \
    --clusters out_dir/TNFa/clusters/chr6_31480875-31598421_cluster.tsv \
    --variant_hashes out_dir/TNFa/variant_hashes.tsv \
    --haploblock_hashes out_dir/haploblock_hashes_chr6.tsv \
    --chr 6 \
    --out out_dir/TNFa/
```

Here is the instruction for **all haploblocks**, please run it separately for every cluster TSV file:
```
python variant_hashes.py \
    --clusters out_dir/clusters/cluster_file.tsv \
    --variant_hashes out_dir/variant_hashes.tsv \
    --haploblock_hashes out_dir/haploblock_hashes_chr6.tsv \
    --chr 6 \
    --out out_dir/
```

This generates variants hashes, 64-character strings of 0/1s. Each individual hash (samples + haplotype) contains:
- strand hash: 4 chars
- chromosome hash: 10 chars
- haploblock hash: 20 chars
- cluster hash: 20 chars
- variant hash: len(variants of interest) chars

The output is a TSV file (with two columns: INDIVIDUAL HASH) is out_dir/individual_hashes.tsv.


## Example Model

Heres an example of how one could take binary strings and feed them into a vector based model with initial weighting...


# Results

## Testing the pipeline

We found 1399 haploblocks in chromosome 6. See [haploblock_boundaries_chr6.tsv](data/haploblock_boundaries_chr6.tsv) for these haploblock boundaries (high recombination rates defined as **rate > 10*average**).

We generated haploblock phased sequences (format: sample_chr_region_start-end_hap0/1.fa) for all samples from the CBH, PUR and GBR populations for the following regions:
- 10 random haploblocks of chr6
- 5 random haploblocks of chr 6
- haploblock overlapping with TNFa
- haploblock overlapping with height genes (TODO ref)

+ haploblock phased sequences and haploblock hashes for TNFa for all populations


# System requirements

During the hackathon we ran the pipeline on a Linux-based machine with 8 CPU cores and 16 GB of memory.


# Acknowlegdements

This work was supported by ELIXIR, the research infrastructure for life science data, and conducted at the ELIXIR BioHackathon Europe.


# References

1. Palsson, G., Hardarson, M.T., Jonsson, H. et al. Complete human recombination maps. Nature 639, 700–707 (2025). https://doi.org/10.1038/s41586-024-08450-5

2. Bjarni V. Halldorsson et al., Characterizing mutagenic effects of recombination through a sequence-level genetic map. Science363,eaau1043 (2019). DOI:10.1126/science.aau1043

3. Yu-Hsiang Tseng, Sumit Walia, Yatish Turakhia, "Ultrafast and ultralarge multiple sequence alignments using TWILIGHT", Bioinformatics, Volume 41, Issue Supplement_1, July 2025, Pages i332–i341, doi: 10.1093/bioinformatics/btaf212
