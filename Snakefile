import os

samples_foldername = '/scratch/users/yzhiyuan/test_snakemake/'
samplenames = ['1001702601_G7']
#input_R1 = ['/scratch/users/yzhiyuan/test_snakemake/1001702601_G7/1001702601_G7_R1_001.fastq.gz',\
#            '/scratch/users/yzhiyuan/test_snakemake/1001702601_G8/1001702601_G8_R1_001.fastq.gz']
#input_R2 = ['/scratch/users/yzhiyuan/test_snakemake/1001702601_G7/1001702601_G7_R2_001.fastq.gz',\
#            '/scratch/users/yzhiyuan/test_snakemake/1001702601_G8/1001702601_G8_R2_001.fastq.gz']
rule all:
  input: expand(samples_foldername + "{sample}/star/Aligned.out.bam", sample=samplenames)
rule STAR:  
  input: reads = ['/scratch/users/yzhiyuan/test_snakemake/1001702601_G7/1001702601_G7_R1_001.fastq.gz',\
                  '/scratch/users/yzhiyuan/test_snakemake/1001702601_G7/1001702601_G7_R2_001.fastq.gz']
  output: "{dir}/{sample}/star/Aligned.out.bam" ,\
          "{dir}/{sample}/star/Unmapped.out.mate1",\
          "{dir}/{sample}/star/Unmapped.out.mate2"
  threads: 6
  params: name="STAR_{sample}", mem="32000"
  log: logout="{dir}/{sample}/log/cluster_{sample}_STAR.log", logerr="{dir}/{sample}/log/cluster_{sample}_STAR.logerr"
  run:
    STAR=os.getenv("STAR", "STAR")
    s=wildcards["sample"]
    pathToGenomeIndex='/scratch/users/yzhiyuan/zhiyuan/ZIKA_pra_Genome2/STAR_DIR'
    check_done=os.path.join(wildcards["dir"],'/', wildcards["sample"], '/', "star.done")
    output_folder=os.path.join(wildcards["dir"], '/', wildcards["sample"], '/', "star/")
    shell("""
        STAR \
        --runThreadN 6 \
        --runMode alignReads \
        --genomeDir {pathToGenomeIndex} \
        --readFilesIn {input.reads[0]} {input.reads[1]} \
        --readFilesCommand zcat \
        --outFilterType BySJout \
        --outFilterMultimapNmax 20 \
        --alignSJoverhangMin 8 \
        --alignSJDBoverhangMin 1 \
        --outFilterMismatchNmax 999 \
        --outFilterMismatchNoverLmax 0.04 \
        --alignIntronMin 20 \
        --alignIntronMax 1000000 \
        --alignMatesGapMax 1000000 \
        --outSAMstrandField intronMotif \
        --outFileNamePrefix {output_folder} \
        --outSAMtype BAM Unsorted \
        --outSAMattributes NH HI AS NM MD \
        --outFilterMatchNminOverLread 0.4 \
        --outFilterScoreMinOverLread 0.4 \
        --clip3pAdapterSeq CTGTCTCTTATACACATCT \
        --outReadsUnmapped Fastx""")
    shell("touch {check_done}")
