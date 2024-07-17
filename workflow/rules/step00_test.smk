rule all:
    input:
        expand("{output_dir}/{sample}.kraken2", output_dir=config["output_dir"], sample=config["samples"])

rule map_to_human:
    input:
        r1="{raw_dir}/{sample}_R1.fastq.gz",
        r2="{raw_dir}/{sample}_R2.fastq.gz",
        ref=config["human_genome"]
    output:
        "{output_dir}/{sample}.sam"
    threads: config["threads"]
    shell:
        """
        bwa mem {input.ref} {input.r1} {input.r2} -M -t {threads} > {output}
        """

rule extract_unmapped:
    input:
        "{output_dir}/{sample}.sam"
    output:
        "{output_dir}/{sample}.unmapped.bam"
    threads: config["threads"]
    shell:
        """
        samtools view -@ {threads} -f 12 -o {output} -bhS {input}
        """

rule sort_bam:
    input:
        "{output_dir}/{sample}.unmapped.bam"
    output:
        "{output_dir}/{sample}.unmapped.sorted.bam"
    threads: config["threads"]
    shell:
        """
        samtools sort -@ {threads} -n {input} -o {output}
        """

rule bam_to_fastq:
    input:
        "{output_dir}/{sample}.unmapped.sorted.bam"
    output:
        r1="{output_dir}/{sample}_unmapped_R1.fastq",
        r2="{output_dir}/{sample}_unmapped_R2.fastq"
    shell:
        """
        bedtools bamtofastq -i {input} -fq {output.r1} -fq2 {output.r2}
        """

rule run_kraken2:
    input:
        r1="{output_dir}/{sample}_unmapped_R1.fastq",
        r2="{output_dir}/{sample}_unmapped_R2.fastq",
        db=config["bacteria_genome"]
    output:
        report="{output_dir}/{sample}.k2report",
        kraken2="{output_dir}/{sample}.kraken2"
    threads: config["threads"]
    shell:
        """
        kraken2 --paired --db {input.db} --threads {threads} \
                --minimum-base-quality 10 --report {output.report} \
                --report-minimizer-data --minimum-hit-groups 3 \
                {input.r1} {input.r2} > {output.kraken2}
        """

