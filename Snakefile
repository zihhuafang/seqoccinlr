import os

include: "rules/common.smk"

SAMPLEFILE = os.path.abspath(config["samples"])
GENOME = os.path.abspath(config["genome"])

samples = get_samples(SAMPLEFILE)
tools = get_tools(config['tools'])

workdir: config['workdir']

wildcard_constraints:
    tool = "(minimap|ngmlr)"

ruleorder: minimap > ngmlr

rule all:
    input:
        expand("stats/reads/{sample}.out", sample=samples.index),
        expand("stats/align/{sample}.{tool}.out", sample=samples.index),
        expand("maping/{sample}.{tool}.bam", sample=samples.index, tool=tools)

##### Modules #####
include: "rules/stats.smk"
include: "rules/mapping.smk"
