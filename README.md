# seqoccinlr
This Snakemake pipeline starts from a set of nanopore reads and an associated genome, simulate reads, map them to the genome with 2 different tools (minimap2 and ngmlr) and compute various statistics about the correctness of the mapping

## Dependencies
  - snakemake-minimal =5.2.4
  - python =3.6.3
  - jinja2 =2.10
  - networkx =2.1
  - matplotlib =2.2.3
  - graphviz =2.38.0
  - bcftools =1.9
  - samtools =1.9
  - bwa =0.7.17
  - pysam =0.15.0
  - nanostat =1.1.0
  - seqtk =1.3
  - minimap2 =2.11
  - nanosim =2.2.0
  - ngmlr =0.2.7
  - drmaa =0.7.6

## Usage

### Step 1: Install workflow

To use this workflow, first download it:
git clone https://github.com/SeqOccin-SV/seqoccinlr.git

This pipeline needs all the tools and versions indicated in the file environment.yaml. An easy way to achieve this is to create a conda environment. For this you need conda (or Miniconda3-4.4.10) and to execute the following commands:
cd seqoccinlr
conda env create --name [yourname] --file environment.yaml

### Step 2: Configure workflow

Configure the workflow according to your needs via editing the file `config.yaml`.

### Step 3: Execute workflow

Test your configuration by performing a dry-run via

    snakemake -n

Execute the workflow locally via

    snakemake --cores $N

using `$N` cores or run it in a cluster environment via

    snakemake --cluster qsub --jobs 100

or

    snakemake --drmaa --jobs 100
