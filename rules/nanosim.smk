
rule subsample:
    input:
        reads=get_fastq
    output:
        "nanosim/{sample}.subset.reads.fa"
    threads:
        1
    params:
        nbsampread=config["nbsampread"]
    log:
        "logs/nanosim/{sample}.subset.log"
    shell:
        "seqtk sample -s$RANDOM {input} {params.nbsampread} | "
        "seqtk seq -a > {output} 2> {log}"

rule read_analysis:
    input:
        "nanosim/{sample}.subset.reads.fa"
    output:
        out1="nanosim/{sample}.training.out",
        out2="nanosim/{sample}.training"
    threads:
        1
    params:
        genome=config["genome"]
    log:
        "logs/nanosim/{sample}.read_analysis.log"
    shell:
        "read_analysis.py -i {input} -r {params.genome} "
        "-o {output.out2} > {output.out1} 2> {log}"

rule unzip:
    input:
        genome=config["genome"]
    output:
        temp("unzip/genome.fa")
    shell:
        "gunzip {input.genome} -c > {output}"

rule simulator:
    input:
        genome="unzip/genome.fa",
        readanal="nanosim/{sample}.training.out",
        train="nanosim/{sample}.training"
    output:
        nanoout="nanosim/{sample}.simulated.out",
        fastaout="nanosim/{sample}.simulated.reads.fa.gz"
    params:
        nbsimread=config["nbsimread"]
    log:
        "logs/nanosim/{sample}.simulator.log"
    shell:
        "simulator.py linear -r {input.genome} -n {params.nbsimread} -c {input.train} "
        "-o {output.nanoout} > {output.nanoout}.log 2> {log};"
        "gzip {output.nanoout}_reads.fasta -c > {output.fastaout}"
