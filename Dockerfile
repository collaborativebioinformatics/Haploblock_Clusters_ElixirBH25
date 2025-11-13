FROM python:3.10-slim

# Set working directory
WORKDIR /app

# Install system dependencies
RUN apt-get update && apt-get install -y --no-install-recommends \
    git \
    wget \
    unzip \
    gzip \
    bzip2 \	
    tar \
    make \
    gcc \
    g++ \
    libbz2-dev \
    libncurses5-dev \
    liblzma-dev \
    libcurl4-openssl-dev \
    && rm -rf /var/lib/apt/lists/*

# -------------------------------------------------------------------
# Install HTSlib, Samtools, and BCFtools from source
# -------------------------------------------------------------------
WORKDIR /opt

# Download source tarballs
RUN wget https://github.com/samtools/samtools/releases/download/1.22.1/samtools-1.22.1.tar.bz2 && \
    wget https://github.com/samtools/bcftools/releases/download/1.22/bcftools-1.22.tar.bz2 && \
    wget https://github.com/samtools/htslib/releases/download/1.22.1/htslib-1.22.1.tar.bz2

# Extract
RUN tar -xvf samtools-1.22.1.tar.bz2 && \
    tar -xvf bcftools-1.22.tar.bz2 && \
    tar -xvf htslib-1.22.1.tar.bz2

# Build and install each tool
RUN cd /opt/htslib-1.22.1 && ./configure && make && make install && \
    cd /opt/samtools-1.22.1 && ./configure && make && make install && \
    cd /opt/bcftools-1.22 && ./configure && make && make install

# Symlinks for global access
RUN ln -s /opt/samtools-1.22.1/samtools /usr/bin/samtools && \
    ln -s /opt/bcftools-1.22/bcftools /usr/bin/bcftools && \
    ln -s /opt/htslib-1.22.1/bgzip /usr/bin/bgzip && \
    ln -s /opt/htslib-1.22.1/tabix /usr/bin/tabix

# -------------------------------------------------------------------
# Install MMSeqs2
# -------------------------------------------------------------------
RUN wget https://mmseqs.com/latest/mmseqs-linux-avx2.tar.gz -O /opt/mmseqs-linux.tar.gz && \
    tar -xvzf /opt/mmseqs-linux.tar.gz -C /opt && \
    ln -s /opt/mmseqs/bin/mmseqs /usr/bin/mmseqs

# -------------------------------------------------------------------
# Clone your repository and set up Python environment
# -------------------------------------------------------------------
WORKDIR /app
#RUN git clone https://github.com/collaborativebioinformatics/Haploblock_Clusters_ElixirBH25.git
RUN git clone https://github.com/MauricioMoldes/Haploblock_Clusters_ElixirBH25.git


# Set repo directory as working directory
WORKDIR /app/Haploblock_Clusters_ElixirBH25

# Copy requirements.txt and install Python dependencies
#COPY requirements.txt .
#RUN pip install --no-cache-dir -r requirements.txt

# -------------------------------------------------------------------
# Download data (optional layer)
# -------------------------------------------------------------------
RUN mkdir -p /app/Haploblock_Clusters_ElixirBH25/data && \
    cd /app/Haploblock_Clusters_ElixirBH25/data && \
    echo "Downloading 1000 Genomes chr6 data..." && \
    wget -q https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000_genomes_project/release/20190312_biallelic_SNV_and_INDEL/ALL.chr6.shapeit2_integrated_snvindels_v2a_27022019.GRCh38.phased.vcf.gz && \
    wget -q https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000_genomes_project/release/20190312_biallelic_SNV_and_INDEL/ALL.chr6.shapeit2_integrated_snvindels_v2a_27022019.GRCh38.phased.vcf.gz.tbi && \
    echo "Downloading chr6 FASTA..." && \
    wget -q https://hgdownload.soe.ucsc.edu/goldenPath/hg38/chromosomes/chr6.fa.gz && \
    gzip -d chr6.fa.gz && \
    bgzip chr6.fa && \
    echo "All data successfully downloaded into /data."

# Default command
CMD ["/bin/bash"]

