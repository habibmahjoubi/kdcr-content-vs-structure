#!/usr/bin/env bash
# Designs 1-3, adventitious-agent panel, reparenting and factorial builds. Idempotent: each stage writes a stamp
# and every output is skipped if present, so the pipeline can be relaunched after a WSL crash.
# Usage: bash pipeline_designs_panel.sh [stage ...]   (default: all stages in order)
set -uo pipefail
K=/home/hm/kdcr
SCR="$(dirname "$(readlink -f "$0")")"; R2="$SCR/../results"
W=$K/work_designs; mkdir -p $W/stamps $W/logs
KR=/home/hm/miniforge3/envs/kdcr/bin/kraken2              # v2.17.1, as in the study
KB=/home/hm/miniforge3/envs/kdcr/bin/kraken2-build
PY=/home/hm/miniconda3/envs/viral-score/bin/python
export PATH=/home/hm/miniforge3/envs/kdcr/bin:$PATH        # kraken2-build helpers (k2mask etc.)
MSL="msl29 msl30 msl31 msl32 msl33 msl34 msl35 msl36 msl37 msl38 msl39 msl40 msl41"
RUNS10="SRR26352207 SRR26352206 SRR26352204 SRR26352205 SRR26352202 SRR26352217 SRR26352216 SRR26352195 SRR26352212 SRR26352208"
NEW4="SRR26352215 SRR26352203 SRR26352184 SRR26352209"   # lab1 RSV 1e4, RSV 1e2, REO 1e4, REO 1e2
RUNS14="$RUNS10 $NEW4"
FQ=$K/data/fastq
stamp() { touch $W/stamps/$1; echo "[$(date +%T)] stage $1 done"; }
done_() { [ -f $W/stamps/$1 ]; }
classify() {  # db out_prefix r1 r2 [extra args]
  local db=$1 out=$2 r1=$3 r2=$4; shift 4
  [ -s ${out}_report.txt ] && return 0
  $KR --db $db --memory-mapping --paired --threads 6 "$@" --report ${out}_report.txt --output ${out}_output.txt $r1 $r2 > /dev/null 2>&1
}

# ---------------------------------------------------------------- S1 missing lab1 runs + pairing check
s1() { done_ s1 && return
  cd $K/data; mkdir -p raw_chunks
  for run in $NEW4; do
    [ -s fastq/${run}_2.fastq ] && continue
    for mate in 1 2; do
      url=$(curl -s "https://www.ebi.ac.uk/ena/portal/api/filereport?accession=${run}&result=read_run&fields=fastq_ftp" | tail -n1 | cut -f2 | tr ';' '\n' | grep "_${mate}.fastq.gz")
      curl -s -r 0-150000000 -o raw_chunks/${run}_${mate}.gz.part "https://${url}"
      { zcat -f raw_chunks/${run}_${mate}.gz.part 2>/dev/null || true; } | head -n 8000000 > fastq/${run}_${mate}.fastq
      rm -f raw_chunks/${run}_${mate}.gz.part
    done
  done
  # mates are synchronised if read i of mate 1 and read i of mate 2 carry the same spot identifier
  { echo -e "run\tpairs_mate1\tpairs_mate2\tidentifier_mismatches"
    for run in $RUNS14; do
      paste <(awk 'NR%4==1{split($1,a,"/"); print a[1]}' fastq/${run}_1.fastq) <(awk 'NR%4==1{split($1,a,"/"); print a[1]}' fastq/${run}_2.fastq) \
      | awk -v r=$run 'BEGIN{n=0;m=0} {n++; if($1!=$2) m++} END{print r"\t"n"\t"n"\t"m}'
    done; } > "$R2/pairing_check.tsv"
  stamp s1; }

# ---------------------------------------------------------------- S2 Design 3 (labels + 13 builds)
s2() { done_ s2 && return
  $PY "$SCR/design3_labels.py" "$R2" > $W/logs/s2_labels.log 2>&1
  for m in $MSL; do
    d=$W/design3/db_$m
    [ -s $d/hash.k2d ] && continue
    mkdir -p $d/taxonomy; cp $K/kraken2_dbs/$m/taxonomy/nodes.dmp $K/kraken2_dbs/$m/taxonomy/names.dmp $d/taxonomy/
    $KB --add-to-library $W/design3/$m.fna --db $d > $W/logs/s2_add_$m.log 2>&1
    $KB --build --db $d --threads 6 > $W/logs/s2_build_$m.log 2>&1
    rm -rf $d/library
  done
  stamp s2; }

# ---------------------------------------------------------------- S3 classify runs (D3 x14, fixed/D1 x new4)
s3() { done_ s3 && return
  for m in $MSL; do
    mkdir -p $W/cls/d3/$m $W/cls/fixed/$m $W/cls/d1/$m
    for r in $RUNS14; do classify $W/design3/db_$m $W/cls/d3/$m/$r $FQ/${r}_1.fastq $FQ/${r}_2.fastq; done
    for r in $NEW4;  do classify $K/kraken2_dbs_fixed/$m $W/cls/fixed/$m/$r $FQ/${r}_1.fastq $FQ/${r}_2.fastq
                        classify $K/kraken2_dbs/$m       $W/cls/d1/$m/$r    $FQ/${r}_1.fastq $FQ/${r}_2.fastq; done
  done
  mkdir -p $W/cls/current
  for r in $NEW4; do classify $K/kraken2_dbs/current $W/cls/current/$r $FQ/${r}_1.fastq $FQ/${r}_2.fastq; done
  stamp s3; }

# ---------------------------------------------------------------- S4 confidence thresholds on classified-read subsets
s4() { done_ s4 && return
  mkdir -p $W/conf/sub
  for r in $RUNS14; do   # reads classified at confidence 0 by the D3 MSL41 database (complete content)
    [ -s $W/conf/sub/${r}_2.fq ] && continue
    awk '$1=="C"{print $2}' $W/cls/d3/msl41/${r}_output.txt > $W/conf/sub/$r.ids
    for mate in 1 2; do
      awk 'NR==FNR{k[$1];next} FNR%4==1{split(substr($1,2),a,"/"); keep=(a[1] in k)} keep' $W/conf/sub/$r.ids $FQ/${r}_${mate}.fastq > $W/conf/sub/${r}_${mate}.fq
    done
  done
  for c in 0.1 0.5; do for design in fixed d1 d3; do for m in $MSL; do
    case $design in fixed) db=$K/kraken2_dbs_fixed/$m;; d1) db=$K/kraken2_dbs/$m;; d3) db=$W/design3/db_$m;; esac
    mkdir -p $W/conf/$design/$c/$m
    for r in $RUNS14; do
      out=$W/conf/$design/$c/$m/$r; [ -s ${out}_report.txt ] && continue
      $KR --db $db --memory-mapping --paired --threads 6 --confidence $c --report ${out}_report.txt --output /dev/null \
          $W/conf/sub/${r}_1.fq $W/conf/sub/${r}_2.fq > /dev/null 2>&1
    done
  done; done; done
  stamp s4; }

# ---------------------------------------------------------------- S5 panel of adventitious agents (simulated)
s5() { done_ s5 && return
  $PY "$SCR/agent_panel.py" simulate > $W/logs/s5_sim.log 2>&1
  P=$W/panel
  for m in $MSL; do
    for design in fixed d1 d3; do
      case $design in fixed) db=$K/kraken2_dbs_fixed/$m;; d1) db=$K/kraken2_dbs/$m;; d3) db=$W/design3/db_$m;; esac
      mkdir -p $P/cls/$design; classify $db $P/cls/$design/$m $P/panel_1.fq $P/panel_2.fq --report-zero-counts
    done
  done
  mkdir -p $P/cls/current; classify $K/kraken2_dbs/current $P/cls/current/current $P/panel_1.fq $P/panel_2.fq --report-zero-counts
  stamp s5; }

# ---------------------------------------------------------------- S6 dated official Kraken2 viral indexes
OFFICIAL="20201202 20210517 20220607 20220908 20221209 20230314 20230605 20231009 20240112 20240605 20240904 20241228 20250402 20250714 20251015 20260226 20260626"
s6() { done_ s6 && return
  O=$W/official; mkdir -p $O/cls
  for d in $OFFICIAL; do
    db=$O/k2_viral_$d
    if [ ! -s $db/hash.k2d ]; then
      mkdir -p $db; curl -s https://genome-idx.s3.amazonaws.com/kraken/k2_viral_$d.tar.gz | tar xz -C $db
    fi
    mkdir -p $O/cls/$d
    for r in $RUNS14; do classify $db $O/cls/$d/$r $FQ/${r}_1.fastq $FQ/${r}_2.fastq; done
    classify $db $O/cls/$d/panel $W/panel/panel_1.fq $W/panel/panel_2.fq --report-zero-counts
  done
  stamp s6; }

# ---------------------------------------------------------------- S7 reparenting, all twelve updates, read-level decomposition
s7() { done_ s7 && return
  $PY "$SCR/reparent_all_updates.py" simulate > $W/logs/s7_sim.log 2>&1
  for m in $MSL; do classify $K/kraken2_dbs_fixed/$m $W/reparent/cls_$m $W/reparent/rep_1.fq $W/reparent/rep_2.fq; done
  stamp s7; }

# ---------------------------------------------------------------- S8 factorial content x taxonomy, and build determinism
s8() { done_ s8 && return
  F=$W/factorial; mkdir -p $F
  # F_a: MSL41 content, current taxonomy.  F_b: MSL41 content + NC_062737.1 (nearest existing ancestor), MSL41 taxonomy.
  # F_c: MSL41 rebuilt identically with 2 threads.  F_d: current rebuilt identically with 2 threads.
  $PY - <<'EOF'
import os
K="/home/hm/kdcr"; F=f"{K}/work_designs/factorial"
src=f"{K}/kraken2_dbs/msl41/viral.filtered.fna"
lab={l.split("\t")[0]:l.split("\t")[2].strip() for l in open(f"{K}/work_designs/design3/msl41_labels.tsv")}
extra="NC_062737.1"
with open(f"{F}/msl41_plus1.fna","w") as fo:
    for l in open(src): fo.write(l)
    keep=False
    for l in open(f"{K}/kraken2_dbs/current/bulk/viral.1.1.genomic.taxid.fna"):
        if l[0]==">":
            keep=l[1:].startswith(extra)
            if keep: fo.write(f">{extra}|kraken:taxid|{lab[extra]} {l.split(' ',1)[1]}")
        elif keep: fo.write(l)
EOF
  build() { local d=$1 fa=$2 tax=$3 th=$4
    [ -s $d/hash.k2d ] && return; mkdir -p $d/taxonomy; cp $tax/nodes.dmp $tax/names.dmp $d/taxonomy/
    $KB --add-to-library $fa --db $d > $d.add.log 2>&1; $KB --build --db $d --threads $th > $d.build.log 2>&1; rm -rf $d/library; }
  build $F/a_msl41content_currenttax $K/kraken2_dbs/msl41/viral.filtered.fna $K/kraken2_dbs/current/taxonomy 6
  build $F/b_msl41plus1_msl41tax     $F/msl41_plus1.fna                       $K/kraken2_dbs/msl41/taxonomy   6
  build $F/c_msl41_rebuild_2threads  $K/kraken2_dbs/msl41/viral.filtered.fna $K/kraken2_dbs/msl41/taxonomy   2
  build $F/d_current_rebuild_2threads $K/kraken2_dbs/current/bulk/viral.1.1.genomic.taxid.fna $K/kraken2_dbs/current/taxonomy 2
  for d in a_msl41content_currenttax b_msl41plus1_msl41tax c_msl41_rebuild_2threads d_current_rebuild_2threads; do
    mkdir -p $F/cls/$d; for r in $RUNS14; do classify $F/$d $F/cls/$d/$r $FQ/${r}_1.fastq $FQ/${r}_2.fastq; done
  done
  { echo -e "database\thash_md5\ttaxid_bits\tcapacity"
    for d in $K/kraken2_dbs/msl41 $K/kraken2_dbs/current $F/a_msl41content_currenttax $F/b_msl41plus1_msl41tax $F/c_msl41_rebuild_2threads $F/d_current_rebuild_2threads; do
      bits=$(grep -ho 'with [0-9]* bits' $d.build.log $K/logs/build_msl41.log 2>/dev/null | head -1 | grep -o '[0-9]*')
      cap=$($KR-inspect --db $d --skip-counts 2>/dev/null | grep 'Table capacity' | grep -o '[0-9]*')
      echo -e "$(basename $d)\t$(md5sum $d/hash.k2d | cut -c1-16)\t$bits\t$cap"
    done; } > "$R2/factorial_databases.tsv"
  stamp s8; }

# ---------------------------------------------------------------- S9 phiX174 reads among classified reads (fixed design)
s9() { done_ s9 && return
  X=$W/phix; mkdir -p $X
  awk '/^>/{p=($0 ~ /^>NC_001422.1/)} p' $K/kraken2_dbs/current/bulk/viral.1.1.genomic.taxid.fna > $X/phix.fna
  for r in SRR26352217 SRR26352207 SRR26352216 SRR26352206 SRR26352195 SRR26352212 SRR26352205 SRR26352208 SRR26352202; do
    [ -s $X/$r.phix_ids ] && continue
    awk '$1=="C"{print $2}' $K/results_fixed/msl30/${r}_output.txt > $X/$r.cls_ids
    awk 'NR==FNR{k[$1];next} FNR%4==1{split(substr($1,2),a,"/"); keep=(a[1] in k); if(keep) print ">"a[1]} FNR%4==2&&keep{print}' $X/$r.cls_ids $FQ/${r}_1.fastq > $X/$r.fa
    blastn -task megablast -query $X/$r.fa -subject $X/phix.fna -evalue 1e-10 -perc_identity 95 -max_hsps 1 -outfmt '6 qseqid' 2>/dev/null | sort -u > $X/$r.phix_ids
  done
  stamp s9; }

# ---------------------------------------------------------------- S10 56-sequence panel: destinations across versions
s10() { done_ s10 && return
  P=$W/p56; mkdir -p $P
  if [ ! -s $P/p56_2.fq ]; then
    for f in $K/positive_control_56/*_1.fastq; do a=$(basename $f _1.fastq); awk -v a=$a 'NR%4==1{sub(/^@/,"@"a"__")}1' $f; done > $P/p56_1.fq
    for f in $K/positive_control_56/*_2.fastq; do a=$(basename $f _2.fastq); awk -v a=$a 'NR%4==1{sub(/^@/,"@"a"__")}1' $f; done > $P/p56_2.fq
  fi
  for m in msl37 msl38 msl39 msl40 msl41; do
    classify $K/kraken2_dbs/$m $P/d1_$m $P/p56_1.fq $P/p56_2.fq
    classify $W/design3/db_$m $P/d3_$m $P/p56_1.fq $P/p56_2.fq
  done
  stamp s10; }

STAGES=${@:-s1 s2 s3 s4 s5 s6 s7 s8 s9 s10}
for s in $STAGES; do echo "[$(date +%T)] start $s"; $s; done
echo "[$(date +%T)] pipeline finished"
