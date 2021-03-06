# Team3-ComparativeGenomics

## Bacterial Comparative Genomics Pipeline

Comparative genomics is a field of biological research in which the genomic features of different organisms are compared. This pipeline is meant to compare different strains that could be clonal and source of the outbreak. Tool and parameter selection is carried out to ensure best performance for de-novo assembled Listeria monocytogenes genomes.

## Installation and Setup 

This pipeline uses as conda based environment to ensure you have the appropriate dependencies. We recommend that you download and install Miniconda from https://conda.io/en/latest/miniconda.html

Example installation for Miniconda on Linux:
```
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
./Miniconda3-latest-Linux-x86_64.sh
rm  Miniconda3-latest-Linux-x86_64.sh
```

Next, clone the repository into your local system:

```
git clone  https://github.gatech.edu/compgenomics2019/Team3-ComparativeGenomics.git
```

Create and activate a conda environment using the yml file provided in our lib folder:

```
#Create environment after downloading yml file
conda-env create -f compgeneff3.yml -n compgeneff3
source activate compgeneff3
```

Install KSNPs using wget, unzip the KSNP package and set KSNP paths:
```
wget https://sourceforge.net/projects/ksnp/files/kSNP3.1_Linux_package.zip
unzip kSNP3.1_Linux_package.zip
vim ~/.bashrc
export PATH=$PATH:$/current_path/kSNP3.1_Linux_package/kSNP3
source ~/.bashrc
vim /current_path/kSNP3.1_Linux_package/kSNP3
set kSNP=/usr/local/kSNP3
set kSNP=/home/lhl/tools/kSNP3.1_Linux_package/kSNP3
export PATH=/current_path/kSNP3.1_Linux_package/kSNP3:$PATH
```

Export path to 'lib' to path variable (installs contains precompiled binaries for pandas, KSNP, Chewbbaca, blast,  which are part of the pipeline)
```
export PATH=$PATH:<path to lib>
```

## Running the pipeline

To run our pipeline with sample data provided in our repository (check sample_input folder)

```
python pipeline.py -m 10 -d genome_aseemblies/ -M
```

For each input genome, the list of generated outputs is as follows:
1. fasta file of assembled genomes. 
