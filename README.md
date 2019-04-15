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
conda-env create -f lib/gp_env.yml -n compgene3
source activate compgene3
```

From within the KSNPs folder in lib, compile binaries for KSNPs (dependency for rnammer which is part of our pipeline):
```
cd KSNPs
./configure
make install
```

Export path to 'lib' to path variable (lib contains precompiled binaries for Chewbbaca, lyve-SET, blast,  which are part of the pipeline)
```
export PATH=$PATH:<path to lib>
```

## Running the pipeline

To run our pipeline with sample data provided in our repository (check sample_input folder)

```
./gp_pipeline.sh -i sample_input -o sample_output
```

For each input genome, the list of generated outputs is as follows:
1. gff file containing the coordinates for the coding sequences
2. fna file for coding nucleotide sequences
3. faa file for coded protein sequences
3. fna file for RNA predictions
