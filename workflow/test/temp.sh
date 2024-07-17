RAW_DIR=
RES_DIR=
THREAD_NUM=
SAMPLE=
H_GENOME_PATH=
M_GENOME_PATH=

###################################################################
#### step 1. human read filtering

# Mapping to Human Genome
bwa mem ${H_GENOME_PATH} \
        ${RAW_DIR}/${SAMPLE}_R1.fastq.gz ${RAW_DIR}/${SAMPLE}_R2.fastq.gz \
        -M  -t ${THREAD_NUM} > ${RES_DIR}/${SAMPLE}.sam

# Extract unmapped paired reads from SAM to BAM
samtools view   -@ ${SAMTOOLS_THREAD} \
                -f 12 -o ${RES_DIR}/${SAMPLE}.unmapped.bam \
                -bhS ${RES_DIR}/${SAMPLE}.sam

# Query sort BAM file
samtools sort -@ ${SAMTOOLS_THREAD} -f -n ${RES_DIR}/${SAMPLE}.unmapped.bam ${RES_DIR}/${SAMPLE}.unmapped.sorted.bam

# Sorted BAM to fastq
bedtools bamtofastq -i ${RES_DIR}/${SAMPLE}.unmapped.sorted.bam \ 
                    -fq ${RES_DIR}/${SAMPLE}_unmapped_R1.fastq -fq2 ${RES_DIR}/${SAMPLE}_unmapped_R2.fastq


###################################################################
#### step 2. Classify microbial reads using kraken2

# Run Kraken2
kraken2 --paired \
        --db ${M_GENOME_PATH} \
        --threads ${THREAD_NUM} \
        --minimum-base-quality 10 \
        --report ${RES_DIR}/${SAMPLE}.k2report \
        --report-minimizer-data \
        --minimum-hit-groups 3 \
        ${RES_DIR}/${SAMPLE}_unmapped_R1.fastq \
        ${RES_DIR}/${SAMPLE}_unmapped_R2.fastq > ${RES_DIR}/${SAMPLE}.kraken2

