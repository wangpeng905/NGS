for sample in samp01 samp02 samp03 samp04 samp05 samp06 samp07 samp08 samp09 samp10
do
# star 1-pass align
mkdir -p /Test170300561/${sample}/Alignment
cd /Test170300561/${sample}/Alignment
/Software/STAR/STAR-2.5.2a/bin/Linux_x86_64/STAR --runThreadN 12 \
	--genomeDir STAR_index \
	--readFilesIn /Test170300561/${sample}/Remove_ribosome_data/${sample}_1.fq.gz /Test170300561/${sample}/Remove_ribosome_data/${sample}_2.fq.gz \
	--readFilesCommand zcat \
	--outFileNamePrefix /Test170300561/${sample}/Alignment/${sample}

# star 2-pass index
mkdir -p /Test170300561/${sample}/STAR2pass_index
cd /Test170300561/${sample}/STAR2pass_index
/Software/STAR/STAR-2.5.2a/bin/Linux_x86_64/STAR --runThreadN 12 \
	--runMode genomeGenerate \
	--genomeDir /Test170300561/${sample}/STAR2pass_index  \
	--genomeFastaFiles hg19.fa \
	--sjdbFileChrStartEnd /Test170300561/${sample}/Alignment/${sample}SJ.out.tab

rm -rf /Test170300561/${sample}/Alignment/${sample}Aligned.out.sam

# star 2-pass align
mkdir -p /Test170300561/${sample}/Alignment_2pass
cd /Test170300561/${sample}/Alignment_2pass
/Software/STAR/STAR-2.5.2a/bin/Linux_x86_64/STAR --runThreadN 12 \
	--genomeDir /Test170300561/${sample}/STAR2pass_index \
	--readFilesIn /Test170300561/${sample}/Remove_ribosome_data/${sample}_1.fq.gz /Test170300561/${sample}/Remove_ribosome_data/${sample}_2.fq.gz \
	--readFilesCommand zcat \
	--outFileNamePrefix /Test170300561/${sample}/Alignment_2pass/${sample}

rm -rf /Test170300561/${sample}/STAR2pass_index/Genome
rm -rf /Test170300561/${sample}/STAR2pass_index/SA
rm -rf /Test170300561/${sample}/STAR2pass_index/SAindex

# picard Add read groups, sort, mark duplicates, and create index
java -jar /Software/picard_tools/picard-tools-2.2.4/picard.jar AddOrReplaceReadGroups \
        I=/Test170300561/${sample}/Alignment_2pass/${sample}Aligned.out.sam \
        O=/Test170300561/${sample}/Alignment_2pass/${sample}_rg_added_sorted.bam \
        SO=coordinate \
        RGID=${sample} \
        RGLB=rna \
        RGPL=illumina \
        RGPU=hiseq \
        RGSM=${sample} 

java -jar /Software/picard_tools/picard-tools-2.2.4/picard.jar MarkDuplicates \
        I=/Test170300561/${sample}/Alignment_2pass/${sample}_rg_added_sorted.bam \
        O=/Test170300561/${sample}/Alignment_2pass/${sample}_dedup.bam  \
        CREATE_INDEX=true \
        VALIDATION_STRINGENCY=SILENT \
        M=/Test170300561/${sample}/Alignment_2pass/${sample}_dedup.metrics

rm -rf /Test170300561/${sample}/Alignment_2pass/${sample}Aligned.out.sam
rm -rf /Test170300561/${sample}/Alignment_2pass/${sample}_rg_added_sorted.bam

# RNA-seq specfic SplitNCigarReads
java -jar GenomeAnalysisTK.jar -T SplitNCigarReads \
        -R hg19.fa \
        -I /Test170300561/${sample}/Alignment_2pass/${sample}_dedup.bam \
        -o /Test170300561/${sample}/Alignment_2pass/${sample}_dedup_split.bam \
        -rf ReassignOneMappingQuality \
        -RMQF 255 \
        -RMQT 60 \
        -U ALLOW_N_CIGAR_READS

rm -rf /Test170300561/${sample}/Alignment_2pass/${sample}_dedup.bam

# indel realign optional
# target_interval file must have file extentsion name likg .bed .list .interval

#java -jar GenomeAnalysisTK.jar \
#    -T RealignerTargetCreator \
#    -R hg19.fa \
#    -known /database/human_db/Mills_and_1000G_gold_standard.indels.hg19.sites.vcf \
#    -I /Test170300561/${sample}/Alignment_2pass/${sample}_dedup_split.bam \
#    -o /Test170300561/${sample}/Alignment_2pass/${sample}_realign_interval.list

#java -jar GenomeAnalysisTK.jar -T IndelRealigner \
#        -R hg19.fa \
#        -I /Test170300561/${sample}/Alignment_2pass/${sample}_dedup_split.bam \
#        -known /database/human_db/Mills_and_1000G_gold_standard.indels.hg19.sites.vcf \
#        -o /Test170300561/${sample}/Alignment_2pass/${sample}_realign.bam \
#        -targetIntervals /Test170300561/${sample}/Alignment_2pass/${sample}_realign_interval.list

# BQSR optional
#java -jar GenomeAnalysisTK.jar \
#        -T BaseRecalibrator \
#        -R hg19.fa \
#        -I /Test170300561/${sample}/Alignment_2pass/${sample}_realign.bam \
#        -knownSites /database/human_db/1000G_phase1.snps.high_confidence.hg19.sites.vcf \
#        -knownSites /database/human_db/Mills_and_1000G_gold_standard.indels.hg19.sites.vcf \
#        -o /Test170300561/${sample}/Alignment_2pass/${sample}_recal_data.table

#java -jar GenomeAnalysisTK.jar  \
#        -T PrintReads \
#        -R hg19.fa  \
#        -I /Test170300561/${sample}/Alignment_2pass/${sample}_realign.bam  \
#        -BQSR /Test170300561/${sample}/Alignment_2pass/${sample}_recal_data.table  \
#        -o /Test170300561/${sample}/Alignment_2pass/${sample}_BQSR.bam 

# variant calling
java -jar GenomeAnalysisTK.jar -T HaplotypeCaller \
        -R hg19.fa \
        -I /Test170300561/${sample}/Alignment_2pass/${sample}_dedup_split.bam \
        -dontUseSoftClippedBases \
        -stand_call_conf 20.0 \
        -o /Test170300561/${sample}/Alignment_2pass/${sample}.vcf

# variant filter
# select SNP 		
java -jar GenomeAnalysisTK.jar \
		-T SelectVariants \
		-R hg19.fa \
		-V /Test170300561/${sample}/Alignment_2pass/${sample}.vcf \
		-selectType SNP \
		-o /Test170300561/${sample}/Alignment_2pass/${sample}.SNP.vcf

# selec INDEL
java -jar GenomeAnalysisTK.jar \
		-T SelectVariants \
		-R hg19.fa \
		-V /Test170300561/${sample}/Alignment_2pass/${sample}.vcf \
		-selectType INDEL \
		-o /Test170300561/${sample}/Alignment_2pass/${sample}.INDEL.vcf

# filter SNP
java -jar GenomeAnalysisTK.jar \
        -T VariantFiltration \
        -R hg19.fa \
        -V /Test170300561/${sample}/Alignment_2pass/${sample}.SNP.vcf \
        -window 35 \
        -cluster 3 \
        -filterName 'MY_SNP_Filter' -filter 'DP < 10.0 || QD < 2.0 || FS > 30.0 || MQ < 40.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0'\
        -o /Test170300561/${sample}/Alignment_2pass/${sample}.SNP.filtered.vcf

rm /Test170300561/${sample}/Alignment_2pass/${sample}.SNP.vc*

# filter INDEL
java -jar GenomeAnalysisTK.jar \
        -T VariantFiltration \
        -R hg19.fa \
        -V /Test170300561/${sample}/Alignment_2pass/${sample}.INDEL.vcf \
        -window 35 \
        -cluster 3 \
        -filterName 'MY_INDEL_Filter' -filter 'DP < 10.0 || QD < 2.0 || FS > 60.0 || ReadPosRankSum < -20.0'\
        -o /Test170300561/${sample}/Alignment_2pass/${sample}.INDEL.filtered.vcf

rm /Test170300561/${sample}/Alignment_2pass/${sample}.INDEL.vc*

# combine SNP and INDEL vcf
java -jar GenomeAnalysisTK.jar \
        -T CombineVariants \
        -R hg19.fa \
        -V /Test170300561/${sample}/Alignment_2pass/${sample}.SNP.filtered.vcf \
        -V /Test170300561/${sample}/Alignment_2pass/${sample}.INDEL.filtered.vcf \
        -o /Test170300561/${sample}/Alignment_2pass/${sample}.ALL.filtered.vcf \
        -genotypeMergeOptions UNSORTED

rm /Test170300561/${sample}/Alignment_2pass/${sample}.SNP.filtered.vc* /Test170300561/${sample}/Alignment_2pass/${sample}.INDEL.filtered.vc*

# filter only passed
vcftools --vcf /Test170300561/${sample}/Alignment_2pass/${sample}.ALL.filtered.vcf --remove-filtered-all --minQ 10 --recode --stdout > /Test170300561/${sample}/Alignment_2pass/${sample}.passed.vcf

rm /Test170300561/${sample}/Alignment_2pass/${sample}.ALL.filtered.vc*

# snpEff annotate
java -jar /home/software/snpEff/snpEff.jar hg19 /Test170300561/${sample}/Alignment_2pass/${sample}.passed.vcf > /Test170300561/${sample}/Alignment_2pass/${sample}.passed.eff.vcf

java -jar /home/software/snpEff/SnpSift.jar annotate -noDownload -v /database/ClinVar/clinvar_20160802+papu.vcf /Test170300561/${sample}/Alignment_2pass/${sample}.passed.eff.vcf > /Test170300561/${sample}/Alignment_2pass/${sample}.passed.eff.clinvar.vcf

java -jar /home/software/snpEff/SnpSift.jar annotate -noDownload -v /database/human_db/dbsnp_138.hg19.vcf /Test170300561/${sample}/Alignment_2pass/${sample}.passed.eff.clinvar.vcf > /Test170300561/${sample}/Alignment_2pass/${sample}.passed.eff.clinvar.dbsnp138.vcf

java -jar /home/software/snpEff/SnpSift.jar annotate -noDownload -v /database/human_db/dbsnp_142.hg19.vcf  /Test170300561/${sample}/Alignment_2pass/${sample}.passed.eff.clinvar.dbsnp138.vcf > /Test170300561/${sample}/Alignment_2pass/${sample}.passed.eff.clinvar.dbsnp138.dbsnp142.vcf 

java -jar /home/software/snpEff/SnpSift.jar annotate -noDownload -v /database/human_db/ESP6500SI-V2-SSA137.GRCh38-liftover.snps_indels.vcf /Test170300561/${sample}/Alignment_2pass/${sample}.passed.eff.clinvar.dbsnp138.dbsnp142.vcf > /Test170300561/${sample}/Alignment_2pass/${sample}.passed.eff.clinvar.dbsnp138.dbsnp142.esp6500.vcf

rm /Test170300561/${sample}/Alignment_2pass/${sample}.passed.eff.vcf
rm /Test170300561/${sample}/Alignment_2pass/${sample}.passed.eff.clinvar.vcf
rm /Test170300561/${sample}/Alignment_2pass/${sample}.passed.eff.clinvar.dbsnp138.vcf
rm /Test170300561/${sample}/Alignment_2pass/${sample}.passed.eff.clinvar.dbsnp138.dbsnp142.vcf
mv /Test170300561/${sample}/Alignment_2pass/${sample}.passed.eff.clinvar.dbsnp138.dbsnp142.esp6500.vcf /Test170300561/${sample}/Alignment_2pass/${sample}.passed.ann.vcf		
done
