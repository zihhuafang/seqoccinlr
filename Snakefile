include: "rules/common.smk"

workdir: config['workdir']

# sample.tsv must be in the working directory
SAMPLEFILE="samples.tsv"

rule all:
    input:
        expand("readstats/{sample}.out", sample=samples.index),
        expand("nanosim/{sample}.simulated.reads.fa.gz", sample=samples.index),
        expand("minimap/{sample}.simulated.reads.genome.aln.sam", sample=samples.index),
        expand("ngmlr/{sample}.simulated.reads.genome.aln.sam", sample=samples.index)


##### Modules #####

include: "rules/readstats.smk"
include: "rules/nanosim.smk"
include: "rules/minimap.smk"
include: "rules/ngmlr.smk"
