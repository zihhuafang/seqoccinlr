
rule read_stats:
    input:
        fastq=get_fastq
    output:
        "readstats/{sample}.out"
    threads:
        4
    resources:
        mem_mb = 8
    log:
        "logs/readstats/{sample}.log"
    shell:
        "NanoStat --fastq {input.fastq} -t {threads} > {output} 2> {log}"
