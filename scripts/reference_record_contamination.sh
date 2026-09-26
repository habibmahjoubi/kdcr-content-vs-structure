#!/usr/bin/env bash
# What are the reads that enter and leave classification between prebuilt databases?
# Samples 200 reads (mate 1) assigned to each high-volume taxon, aligns them to the RefSeq viral library
# (megablast, subject coordinates) and to the Chinese hamster genomes (minimap2).
set -uo pipefail
K=/home/hm/kdcr; W=$K/work_exposure/contam; mkdir -p $W
E=/home/hm/miniforge3/envs/kdcr/bin; export PATH=$E:$PATH
FQ=$K/data/fastq
pick() {  # output_file taxid tag run
  awk -v t=$2 '$1=="C" && $3==t {print $2}' $1 | shuf -n 200 --random-source=<(yes) > $W/$3_$4.ids
  awk 'NR==FNR{k[$1];next} FNR%4==1{split(substr($1,2),a,"/"); keep=(a[1] in k); if(keep) print ">"a[1]} FNR%4==2&&keep{print}' \
      $W/$3_$4.ids $FQ/${4}_1.fastq > $W/$3_$4.fa
}
for run in SRR26352206 SRR26352207; do
  pick $K/work_designs/official/cls/20220607/${run}_output.txt 1969841 proteus_isfahan $run
  pick $K/work_designs/official/cls/20220607/${run}_output.txt 2500149 salmonella_TS13 $run
  pick $K/work_designs/official/cls/20220607/${run}_output.txt 159150  shamonda $run
done
cat $W/*.fa > $W/all.fa
blastn -task megablast -query $W/all.fa -db $K/work/blastdb/viral_current -evalue 1e-5 -max_target_seqs 1 -max_hsps 1 \
       -outfmt '6 qseqid sseqid pident length sstart send evalue' -num_threads 6 > $W/all_viral.tsv 2>/dev/null
for g in picr chok1; do
  minimap2 -x sr -c --secondary=no -I 1G -t 6 $K/work_exposure/host/$g.fna $W/all.fa > $W/all_$g.paf 2>/dev/null
done
echo done
