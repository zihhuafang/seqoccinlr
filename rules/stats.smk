
rule read_stats:
    input:
        fastq=get_fastq
    output:
        "stats/erads/{sample}.out"
    threads:
        4
    resources:
        mem_mb = 16
    log:
        "logs/stats/reads/{sample}.{tool}.log"
    shell:
        "NanoStat --fastq {input.fastq} -t {threads} > {output} 2> {log}"

rule align_stats:
    input:
        bam="mapping/{sample}.{tool}.bam"
    output:
        "stats/align/{sample}.{tool}.out"
    threads:
        4
    resources:
        mem_mb = 16
    log:
        "logs/stats/align/{sample}.{tool}.log"
    shell:
        "NanoStat --bam {input.bam} -t {threads} > {output} 2> {log}"
