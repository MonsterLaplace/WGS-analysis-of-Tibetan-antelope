#!/usr/bin/bash

## 这是多样本变异检测流程的后半部分，这个部分假设你已经用fastq_to_gvcf.sh为每个样本独立生成了其对应的比对结果和gvcf

#reference
reference=/mnt/z/xb/TAntelope/00.rawdata/genome/TAntelope.fa

## shell执行参数
samples=$1  ## 所有的样本ID，用","分开
indir=$2  ## 输入目录的路径，这个输入路径要与fastq_to_gvcf.sh的输出路径完全相同
outname=$3 ## 设置输出文件名的前缀

outdir=$indir ## 输入和输出路径相同

# 设置群体变异检测结果的输出目录
if [ ! -d $outdir/population ]
then mkdir -p $outdir/population
fi

## 按照","，把所有样本ID拆分出来存至数组中
samples=$(echo $samples | tr "," "\n")

# ## 第二，按照染色体分开，合并每个样本对应染色体的gVCF，然后再各自对染色体进行Joint calling(注意，这里和GATK3不同，GATK4只能接受一个gvcf的输入，因此需要先合并),
# ## 最后合并各个染色体的Genotype结果, 这要求每个样本必须按照染色体输出gvcf，速度快
 chrom=( Chr01 Chr02 Chr03 Chr04 Chr05 Chr06 Chr07 Chr08 Chr09 Chr10 Chr11 Chr12 Chr13 Chr14 Chr15 Chr16 Chr17 Chr18 Chr19 Chr20 Chr21 Chr22 Chr23 Chr24 Chr25 Chr26 Chr27 Chr28 Chr29 ChrX ChrY )

 for i in ${chrom[@]}; do
     sample_gvcfs=""
     for sample in $samples ; do
 	    sample_gvcfs=${sample_gvcfs}"-V $outdir/${sample}/gatk/${sample}.HC.${i}.g.vcf.gz  "
     done
     time gatk CombineGVCFs \
         -R $reference \
         ${sample_gvcfs} \
         -O $outdir/population/${outname}.HC.${i}.g.vcf.gz && echo "** ${outname}.HC.${i}.g.vcf.gz done ** " && \
     time gatk GenotypeGVCFs \
         -R $reference \
         -V $outdir/population/${outname}.HC.${i}.g.vcf.gz \
         -O $outdir/population/${outname}.HC.${i}.vcf.gz && echo "** ${outname}.HC.${i}.vcf.gz done ** " &
 done && wait
 merge_vcfs=""
 for i in ${chrom[@]}; do
     merge_vcfs=${merge_vcfs}"-I $outdir/population/${outname}.HC.${i}.vcf.gz  "
 done && time gatk MergeVcfs ${merge_vcfs} -O $outdir/population/${outname}.merged.HC.vcf.gz && echo "** MergeVcfs done **"

# 使用SelectVariants，选出SNP,为SNP作硬过滤
time gatk SelectVariants \
     -select-type SNP \
	 -V $outdir/population/${outname}.merged.HC.vcf.gz \
	 -O $outdir/population/${outname}.HC.snp.vcf.gz && \
time gatk VariantFiltration \
     -V $outdir/population/${outname}.HC.snp.vcf.gz \
	 --filter-expression "QD < 2.0 || MQ < 40.0 || FS > 60.0 || SOR > 3.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0" \
	 --filter-name "Filter" \
	 -O $outdir/population/${outname}.snp.filter.vcf.gz && echo "** SNPs filter done **"

# 使用SelectVariants，选出Indel,为Indel作硬过滤
time gatk SelectVariants \
     -select-type INDEL \
	 -V $outdir/population/${outname}.merged.HC.vcf.gz \
	 -O $outdir/population/${outname}.HC.indel.vcf.gz && \
time gatk VariantFiltration \
     -V $outdir/population/${outname}.HC.indel.vcf.gz \
	 --filter-expression "QD < 2.0 || FS > 200.0 || SOR > 10.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0" \
	 --filter-name "Filter"  \
	 -O $outdir/population/${outname}.indel.filter.vcf.gz && echo "** indels filter done **"

# 重新合并过滤后的SNP和Indel
#time gatk MergeVcfs \
     -I $outdir/population/${outname}.snp.filter.vcf.gz \
	 -I $outdir/population/${outname}.indel.filter.vcf.gz \
	 -O $outdir/population/${outname}.filter.vcf.gz && echo "** SNPs and Indels filter done **"

# 可以被删除清理的文件，这不是必须执行的,可以保留.g.vcf.gz，原始HC.vcf.gz和完成质控的HC.filter.vcf.gz
 rm -f $outdir/population/${outname}.HC.vcf.gz 
 $outdir/population/${outname}.HC.*.recal