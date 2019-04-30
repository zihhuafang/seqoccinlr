import os

include: "rules/common.smk"

SAMPLEFILE = os.path.abspath(config["samples"])
GENOME = os.path.abspath(config["genome"])

samples = get_samples(SAMPLEFILE)
tools = get_tools(config['tools'])

workdir: config['workdir']
wildcard_constraints:
    tool = "(minimap|ngmlr)"

rule all:
    input:
        expand("readstats/{sample}.out", sample=samples.index),


##### Modules #####

include: "rules/stats.smk"
include: "rules/mapping.smk"
