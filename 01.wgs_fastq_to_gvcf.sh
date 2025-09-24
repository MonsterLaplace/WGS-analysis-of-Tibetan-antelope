#!/usr/bin/bash

## 这是多样本变异检测流程的前半部分，这个部分假设你的每个样本只有一对用Illumina测序仪测序的PE fastq

#reference
reference=/mnt/z/xb/TAntelope/00.rawdata/genome/TAntelope.fa

## shell执行参数
fq1=$1
fq2=$2
RGID=$3  ## Read Group，一般用Lane ID代替
library=$4  ## 测序文库编号
sample=$5  ## 样本ID
outdir=$6  ## 输出目录的路径

## 按样本设置目录
outdir=${outdir}/${sample}

## 通过fastq1获得fastq的前缀名字，这里假设了原始的fastq1和fastq2有相同的前缀名字
## 并且假设fastq1的文件名格式为*.1.fq.gz;
file=`basename $fq1`
base=${file%.1.fq.gz}

# output diretory
if [ ! -d $outdir/bwa ]
then mkdir -p $outdir/bwa
fi

if [ ! -d $outdir/gatk ]
then mkdir -p $outdir/gatk
fi

if [ ! -d $outdir/cleanfq ]
then mkdir -p $outdir/cleanfq
fi

if [ ! -d $outdir/cleanfq/report ]
then mkdir -p $outdir/cleanfq/report
fi

#质控
fastp \
-i $fq1 \
-I $fq2 \
-o $outdir/cleanfq/${base}.cleanR1.fq.gz \
-O $outdir/cleanfq/${base}.cleanR2.fq.gz \
-W 5 \
-M 20 \
-5 \
-q 15 \
-u 40 \
-n 0 \
-l 75 \
-w 5 \
-j $outdir/cleanfq/report/${fq_file_name}.clean.json \
-h $outdir/cleanfq/report/${fq_file_name}.clean.html \
-R $outdir/cleanfq/report/${fq_file_name}.report && echo "** fq QC done **"
  
##转录组使用bwa比对并排序
bwa-mem2 mem -t 250 -R '@RG\tID:$RGID\tPL:ILLUMINA\tLB:$library\tSM:$sample' $reference \
	$outdir/cleanfq/${base}.cleanR1.fq.gz $outdir/cleanfq/${base}.cleanR2.fq.gz | samtools view -@ 60 -m 32G -Sb - > $outdir/bwa/${sample}.bam && \
	echo "** BWA MEM done **"

samtools sort -@ 8 -m 32G -O bam -o $outdir/bwa/${sample}.sorted.bam $outdir/bwa/${sample}.bam && echo "** sorted raw bamfile done "

##修改BAM文件的Read Group
gatk AddOrReplaceReadGroups \
  -I $outdir/bwa/${sample}.sorted.bam \
  -O $outdir/bwa/${sample}.addgroup.sorted.bam \
  --RGID $RGID \
  --RGLB $library \
  --RGPL ILLUMINA \
  --RGPU unit1 \
  --RGSM $sample \
  --VALIDATION_STRINGENCY LENIENT && \
  echo "** add read group to bamfile done **"

## 这一步不是必须的 
samtools index $outdir/bwa/${sample}.addgroup.sorted.bam && echo "** ${sample}.sorted.bam index done **"

## 标记重复序列 
gatk MarkDuplicates \
  --TMP_DIR /mnt/z/xb/TAntelope/03.cache \
  -I $outdir/bwa/${sample}.addgroup.sorted.bam \
  -O $outdir/bwa/${sample}.sorted.markdup.bam \
  -M $outdir/bwa/${sample}.markdup_metrics.txt && echo "** ${sample}.sorted.bam MarkDuplicates done **"
  
## 为${sample}.sorted.markdup.bam构建Index，这是继续后续步骤所必须的
samtools index $outdir/bwa/${sample}.sorted.markdup.bam && echo "** ${sample}.sorted.markdup.bam index done **"

## 分染色体call每个样品的SNP
chrom=( Chr01 Chr02 Chr03 Chr04 Chr05 Chr06 Chr07 Chr08 Chr09 Chr10 Chr11 Chr12 Chr13 Chr14 Chr15 Chr16 Chr17 Chr18 Chr19 Chr20 Chr21 Chr22 Chr23 Chr24 Chr25 Chr26 Chr27 Chr28 Chr29 ChrX ChrY )
 for i in ${chrom[@]};do
 	time gatk HaplotypeCaller \
      --native-pair-hmm-threads 10 \
      --emit-ref-confidence GVCF \
	 -R $reference \
	 -I $outdir/bwa/${sample}.sorted.markdup.bam -L $i \
 	 -O $outdir/gatk/${sample}.HC.$i.g.vcf.gz && echo "** ${sample}.HC.$i.g.vcf.gz done **" &
 done

 ## 删除过程文件
rm -f $outdir/bwa/${sample}.bam $outdir/bwa/${sample}.sorted.bam $outdir/bwa/${sample}.sorted.bam.bai
