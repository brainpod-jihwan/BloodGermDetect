import os

# Config
configfile: "config.yaml"

# Command line parameters
output_dir = config["output_dir"]
output_prefix = config["output_prefix"]
SAMPLE = config["sample"]
UCSC_GENOME = config["ucsc_genome"]
RAW_DIR = config["raw_dir"]
THREAD_NUM = config["thread_num"]
SAMTOOLS_THREAD = config["samtools_thread"]


#####################################
# Rule to ensure the output directory exists
rule create_output_dir:
    output:
        directory(output_dir)
    shell:
        """
        mkdir -p {output}
        """

#####################################
# Rule to map reads to Human Genome
rule map_to_human_genome:
    input:
        R1 = f"{RAW_DIR}/{SAMPLE}_R1.fastq.gz",
        R2 = f"{RAW_DIR}/{SAMPLE}_R2.fastq.gz"
    output:
        sam = os.path.join(output_dir, f"{output_prefix}.F1.sam")
    log:
        os.path.join(output_dir, "01.Align_to_UCSC_genome.log")
    threads: THREAD_NUM
    shell:
        """
        bwa mem {UCSC_GENOME} {input.R1} {input.R2} -M -t {threads} > {output.sam} 2> {log}
        """

#####################################
# Rule to extract unmapped reads
rule extract_unmapped_reads:
    input:
        sam = os.path.join(output_dir, f"{output_prefix}.F1.sam")
    output:
        bam = os.path.join(output_dir, f"{output_prefix}.unmapped1.bam"),
        sorted_bam = os.path.join(output_dir, f"{output_prefix}.unmapped1.sorted.bam"),
        R1_fastq = os.path.join(output_dir, f"{output_prefix}_F1_R1.fastq"),
        R2_fastq = os.path.join(output_dir, f"{output_prefix}_F1_R2.fastq")
    log:
        os.path.join(output_dir, "02.filter_unmapped_reads1.log")
    threads: SAMTOOLS_THREAD
    shell:
        """
        samtools view -@ {threads} -f 4 -F 256 -o {output.bam} -bhS {input.sam} 2> {log}
        samtools sort -@ {threads} -n {output.bam} -o {output.sorted_bam} 2>> {log}
        bedtools bamtofastq -i {output.sorted_bam} -fq {output.R1_fastq} -fq2 {output.R2_fastq} 2>> {log}
        rm {output.bam} {input.sam}
        """

#####################################
# Workflow
rule all:
    input:
        os.path.join(output_dir, f"{output_prefix}_F1_R1.fastq"),
        os.path.join(output_dir, f"{output_prefix}_F1_R2.fastq")
