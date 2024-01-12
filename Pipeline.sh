## Title     : PIPELINE
## Objective : Workflow for virus integration detection
## Created by: fpeng
## Created on: 2024/01/12



## extract barcode from BD sequencing ##
## used with umi_tools ##
/public/home/fpeng/miniconda3/envs/scrna/bin/umi_tools whitelist --stdin /public/home/fpeng/single_cell_rawdata/con1_R1.fastq.gz \
        --extract-method=regex \
        --bc-pattern="(?P<cell_1>.{9})(?P<discard_1>.{12})(?P<cell_2>.{9})(?P<discard_2>.{13})(?P<cell_3>.{9})(?P<umi_1>.{8})T+" \
        --plot-prefix=/public/home/fpeng/single_cell_rawdata/testfile/con1_true2 --log2stderr --method=umis \
        --knee-method=density --allow-threshold-error > /public/home/fpeng/single_cell_rawdata/testfile/con1_BDwhitelist.txt

/public/home/fpeng/miniconda3/envs/scrna/bin/umi_tools extract --bc-pattern="(?P<cell_1>.{9})(?P<discard_1>.{12})(?P<cell_2>.{9})(?P<discard_2>.{13})(?P<cell_3>.{9})(?P<umi_1>.{8})T+" --stdin /public/home/fpeng/single_cell_rawdata/con1_R1.fastq.gz --stdout /public/home/fpeng/single_cell_rawdata/testfile/con1_R1_extracted.fastq.gz --extract-method=regex --read2-in /public/home/fpeng/single_cell_rawdata/con1_R2.fastq.gz --read2-out=/public/home/fpeng/single_cell_rawdata/testfile/con1_R2_extracted.fastq.gz --error-correct-cell --whitelist=/public/home/fpeng/single_cell_rawdata/testfile/con1_BDwhitelist.txt


## align with genomes ##
## used with STAR ##
STAR --runMode genomeGenerate --genomeDir STAR_index --runThreadN 5 --genomeFastaFiles Gh38_hpvall.fa --sjdbGTFfile Gh38_hpvall.gtf --sjdbOverhang 149
/public/home/fpeng/miniconda3/envs/scrna/bin/STAR --runThreadN 5 --genomeDir /public/home/fpeng/ncbi/public/sra/STAR_HPV_GH38/STAR_index --readFilesIn /public/home/fpeng/single_cell_rawdata/testfile/con5_R2_extracted.fastq.gz --rand zcat --chimOutType Junctions SeparateSAMold WithinBAM HardClip --chimSegmentMin 20 --chimJunctionOverhangMin 20 --outSAMtypByCoordinate --outFileNamePrefix /public/home/fpeng/ncbi/public/sra/STAR_HPV_GH38/con5_star/con5_star --outTmpDir /public/home/ublic/sra/STAR_HPV_GH38/con5_star/Temp


## extract id from juction result ##
## used shell ##
awk '{print $10}' case3.junction > case3_junction.ids

## from sam to bam ##
## used samtools ##
../../../seqtools/samtools-1.14/samtools view -@6 -bS case3_starChimeric.out.sam > case3_starChimeric.out.bam

## extract bam from chimeric result ##
## used picard ##
java -jar ../../../seqtools/picard.jar FilterSamReads I=con5_starChimeric.out.bam O=con5_Chimeric.bam READ_LIST_FILE=con5_junction.id FILTER=includeReadList

## from bam to fastq ##
## used samtools ##
../../../seqtools/samtools-1.14/samtools fastq con5_Chimeric.bam > con5_Chimeric.fq

## from fastq to fasta ##
## used shell ##
awk ' {if (NR%4 == 1) {print ">" substr ($0, 2)}} {if (NR%4 == 2) {print}}' con5_Chimeric.fq > con5_Chimeric.fa


## make blast db ##
## used blast ##
makeblastdb -in ../../GRCh38.p13.genome.fa -dbtype nucl -parse_seqids -out ./GH38
makeblastdb -in ../../hpvall.fa -dbtype nucl -parse_seqids -out ./hpvall


## blast align ##
## used  blast ##
blastn -query con5_Chimeric.fa -db ../blast_index/GH38 -evalue 1e-10 -outfmt 6 -num_threads 6 -out con5_GH38.blast
blastn -query con5_Chimeric.fa -db ../blast_index/hpvall -evalue 1e-5 -outfmt 6 -num_threads 6 -out con5_hpvall.blast


## merge blast result ##
## used shell ##
sort -t $'\t' -k 1,1 -u case3_GH38.blast > case3_GH38.sort.blast
sort -t $'\t' -k 1,1 -u case3_hpvall.blast > case3_hpvall.sort.blast
cat case3_hpvall.sort.blast case3_GH38.sort.blast | sort | uniq -w 75 -D > case2_GH38_hpvall.blast
sort -t $'\t' -k 1,1 -k 7n,7 case12_GH38_hpvall.blast > case12_GH38_hpvall.sort.blast

## find breakpoints ##
## used sites_cluster.py ##
python ../sites_cluster.py case11_GH38_hpvall.sort.blast

## remove duplicated ##
## used shell and only_bk.py ##
awk -v OFS="\t" '{if(NR>1) {print $2,$3,$4,$5}}' sites_cluster_result| sort -t $'\t' -k 1,1 -k 2,2n > only_bk.txt
python ../only_bk.py only_bk.txt > only_bk_result.txt



