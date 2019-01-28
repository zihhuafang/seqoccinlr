
rule minimap:
    input:
        "nanosim/{sample}.simulated.reads.fa.gz"
    output:
        "minimap/{sample}.simulated.reads.genome.aln.sam"
    threads:
        4
    resources:
        mem_mb=30
    params:
        genome=config["genome"]
    log:
        "logs/minimap/{sample}.simulated.reads.genome.aln.log"
    shell:
        "minimap2 -t {threads} -ax map-ont {params.genome} {input} > {output} 2> {log}"
