Chip-seq
# cd /mnt/Storage/home/xiuwenchao/MTS/new_chip-seq/$1
# bowtie2 -x /mnt/Storage/home/xiuwenchao/bin/bowtie2/mm9 -1 ./$1.1.fq -2 ./$1.2.fq -S $1.sam
# samtools view -bht /mnt/Storage/home/xiuwenchao/bin/chromInfo_mm9.txt $1.sam > $1.bam
# /mnt/Storage/home/wangcf/bin/bamtools/bin/bamtools filter -mapQuality ">=30" -in $1.bam -out $1.filtered.bam
# cd /mnt/Storage/home/xiuwenchao/MTS/TS_oeEsrrb_Chip-seq/$1
# macs2 callpeak -g mm -n $1 -B -q 0.05 --nomodel --shift=73 --SPMR -t $1.filtered.bam -c /mnt/Storage/home/xiuwenchao/MTS/TS_oeEsrrb_Chip-seq/Input_2nd/Input_2nd.filtered.bam
# cd /mnt/Storage/home/xiuwenchao/MTS/new_chip-seq/$1
# bedtools intersect -a $1_treat_pileup.bdg -b /mnt/Storage/home/xiuwenchao/bin/chr_limit_mm9.bed -wa -f 1.00 > $1_treat_pileup.bdg.tmp
# bedGraphToBigWig $1_treat_pileup.bdg.tmp /mnt/Storage/home/xiuwenchao/bin/mm9.len $1.bw

RNA-seq
# cd /mnt/Storage/home/xiuwenchao/MTS/new_RNA_seq2
# tophat -p 8 --mate-inner-dist 0 --mate-std-dev 20 -o ./$1 -G /mnt/Storage/home/huse/Data/RNA_pipe/mm9_refSeq.gtf /mnt/Storage/data/Bowtie/mm9 ./$1/$1.1.fastq ./$1/$1.2.fastq &
# read_distribution.py -i ./accepted_hits.bam -r /mnt/Storage/home/xiuwenchao/bin/mm9.refseq.bed
# htseq-count -f bam ./accepted_hits.bam -i transcript_id /mnt/Storage/home/wangcf/annotations/mm9_refSeq.gtf > $1.count
# cufflinks -G /mnt/Storage/home/wangcf/annotations/mm9_refSeq.gtf --output-dir ./ ./accepted_hits.bam

bw correlation
bigwig_correlation.py -r EGFP_newmethod.r -f PDF TSC_EGFP.bw CJY_1008_13.bw &

'UCSC track'
track type=bigWig name="WT_mc" description="WT_mc" bigDataUrl=http://compbio.tongji.edu.cn/~xiuwenchao/WT_mc.bw
track type=bigWig name="WT_hmc" description="WT_hmc" bigDataUrl=http://compbio.tongji.edu.cn/~xiuwenchao/WT_hmc.bw
track type=bigWig name="Tet2_homo_mc_1" description="Tet2_homo_mc_1" bigDataUrl=http://compbio.tongji.edu.cn/~xiuwenchao/Tet2_homo_mc_1.bw
track type=bigWig name="Tet2_homo_mc_2" description="Tet2_homo_mc_2" bigDataUrl=http://compbio.tongji.edu.cn/~xiuwenchao/Tet2_homo_mc_2.bw
track type=bigWig name="Tet2_homo_hmc_1" description="Tet2_homo_hmc_1" bigDataUrl=http://compbio.tongji.edu.cn/~xiuwenchao/Tet2_homo_hmc_1.bw
track type=bigWig name="Tet2_homo_hmc_2" description="Tet2_homo_hmc_2" bigDataUrl=http://compbio.tongji.edu.cn/~xiuwenchao/Tet2_homo_hmc_2.bw
