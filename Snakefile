configfile: "config.yaml"

rule all:
    input:
        "read_stats/nanostat.out",
        "read_sample/read.sample.fastq",
        "data/reads/reads.fa",
        "nanosim/read_analysis.out"

rule read_stats:
    input:
        fa="data/genome/genome.fa.gz",
        fq="data/reads/reads.fastq.gz"
    output:
        "read_stats/nanostat.out"
    threads:
        4
    resources:
        mem_mb=8
    log:
        "logs/read_stats/nanostat.log"
    shell:
        "NanoStat --fastq {input.fq} -t {threads} > {output} 2> {log}"

rule read_sample:
    input:
        "data/reads/reads.fastq.gz"
    output:
        "read_sample/read.sample.fastq"
    threads:
        1
    params:
        nbreads=config["nbreads"]
    log:
        "logs/read_sample/seqtk.sample.log"
    shell:
        "seqtk sample -s$RANDOM {input} {params.nbreads} > {output} 2> {log}"

rule convert:
    input:
        "data/reads/reads.fastq.gz"
    output:
        "data/reads/reads.fa"
    threads:
        1
    log:
        "logs/convert/seqtk.convert.log"
    shell:
        "seqtk seq -a {input} > {output} 2> {log}"

rule read_analysis:
    input:
        reads="data/reads/reads.fa",
        genome="data/genome/genome.fa.gz"
    output:
        "nanosim/read_analysis.out"
    threads:
        1
    log:
        "logs/nanosim/read_analysis.log"
    shell:
        "read_analysis.py -i {input.reads} -r {input.genome} > {output} 2> {log}"
