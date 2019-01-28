
rule ngmlr:
    input:
        "nanosim/{sample}.simulated.reads.fa.gz"
    output:
        "ngmlr/{sample}.simulated.reads.genome.aln.sam"
    threads:
        4
    resources:
        mem_mb=30
    params:
        genome=config["genome"]
    log:
        "logs/ngmlr/{sample}.simulated.reads.genome.aln.log"
    shell:
        "ngmlr -t {threads} -r {params.genome} -q {input} -o {output} -x ont 2> {log}"
