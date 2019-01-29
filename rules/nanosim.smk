
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
        out2="nanosim/{sample}.training_model_profile"
    threads:
        1
    params:
        genome=config["genome"],
        prefix="nanosim/{sample}.training"
    log:
        "logs/nanosim/{sample}.read_analysis.log"
    shell:
        "read_analysis.py -i {input} -r {params.genome} "
        "-o {params.prefix} > {output.out1} 2> {log}"

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
        train="nanosim/{sample}.training_model_profile"
    output:
        nanoout="nanosim/{sample}.simulated.out",
        fastaout="nanosim/{sample}.simulated.reads.fa.gz"
    params:
        nbsimread=config["nbsimread"],
        prefix1="nanosim/{sample}.training",
        prefix2="nanosim/{sample}.simulating"
        #lambda w: "nanosim/%s.training" %(w.sample)
    log:
        "logs/nanosim/{sample}.simulator.log"
    shell:
        "simulator.py linear -r {input.genome} -n {params.nbsimread} -c {params.prefix1} "
        "-o {params.prefix2} > {output.nanoout} 2> {log};"
        "gzip {params.prefix2}_reads.fasta -c > {output.fastaout}"
