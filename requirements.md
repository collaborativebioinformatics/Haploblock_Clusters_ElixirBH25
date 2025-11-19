### Install Python libraries via [Python venv](https://docs.python.org/3/library/venv.html) with the following command:

```
python3 -m venv --system-site-packages ~/pyEnv_ElixirBH2025
source ~/pyEnv_ElixirBH2025/bin/activate
pip install --upgrade pip
pip install numpy pyyaml
```

### Install samtools, tabix and bcftools
Download software from https://www.htslib.org/
```
wget https://github.com/samtools/samtools/releases/download/1.22.1/samtools-1.22.1.tar.bz2
wget https://github.com/samtools/bcftools/releases/download/1.22/bcftools-1.22.tar.bz2
wget https://github.com/samtools/htslib/releases/download/1.22.1/htslib-1.22.1.tar.bz2
```

Note: samtools requires also
```
sudo dnf install ncurses-devel (Redhat)
```
or
```
sudo apt install libncurses5-dev liblzma-dev libcurl4-openssl-dev (Ubuntu or Debian)
```

Note: bcftools requires also
```
sudo apt install libbz2-dev (Ubuntu)
```

Then do:
```
tar -xvf samtools-1.22.1.tar.bz2
tar -xvf bcftools-1.22.tar.bz2
tar -xvf htslib-1.22.1.tar.bz2
```

```
cd samtools-1.22.1/
./configure
make
make install
```

```
cd bcftools-1.22/
./configure
make
make install
```

```
cd htslib-1.22.1/
./configure
make
make install
```

###  Create simlinks for access everywhere (or export PATH)
```
cd /usr/bin
sudo ln -s ~/samtools-1.22.1/samtools .
sudo ln -s ~/bcftools-1.22/bcftools .
sudo ln -s ~/htslib-1.22.1/bgzip .
sudo ln -s ~/htslib-1.22.1/tabix .
```

### Install MMSeqs2 (https://github.com/soedinglab/MMseqs2.git)
Download a binary (CPU or GPU) from https://dev.mmseqs.com/latest/
```
tar -xvzf mmseqs-linux-*.tar.gz
cd /usr/bin
ln -s mmseqs/bin/mmseqs
```
