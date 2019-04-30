

rule ngmlr:
    input:
        reads=get_fastq
    output:
        bam="mapping/{sample}.ngmlr.bam"
    threads:
        20
    resources:
        mem_mb=30
    params:
        genome=GENOME
    log:
        "logs/mapping/{sample}.ngmlr.log"
    shell:
        "ngmlr -t {threads} -r {params.genome} -q {input.reads} "
        "-x ont | "
        "samtools view -bS | "
        "samtools sort -@{threads} -o {output.bam} 2> {log}"

rule minimap:
    input:
        reads=get_fastq
    output:
        bam="mapping/{sample}.minimap.bam"
    threads:
        20
    resources:
        mem_mb=30
    params:
        genome=GENOME
    log:
        "logs/minimap/{sample}.log"
    shell:
        "minimap2 -t {threads} -ax map-ont {params.genome} {input} |
        "samtools view -bS | "
        "samtools sort -@{threads} -o {output.bam} 2> {log}"


rule samtobam:
    input:
        sam="mapping/{sample}.{tool}.sam"
    output:
        bam="mapping/{sample}.{tool}.bam"
    threads:
        4
    log:
        "logs/samtobam/{sample}.{tool}.log"
    shell:
        "samtools view -bS {input.sam} | "
        "samtools sort -@{threads} -o {output.bam} 2> {log}"
