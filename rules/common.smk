import pandas as pd


###### Config file and sample sheets #####
configfile: "config.yaml"
workdir: config['workdir']

SAMPLEFILE="samples.tsv"
samples = pd.read_table(SAMPLEFILE).set_index("sample", drop=False)

def get_fastq(wildcards):
    """Get fastq files of given sample-unit."""
    fastq = samples.loc[(wildcards.sample), ["fastq"]].dropna()
    return fastq
