# kdcr-content-vs-structure

Code, scripts, and summary result tables for:

> Mahjoubi H. **Structural Change Is Not Impact: Why Reference-Artifact Diffs Fail to
> Scope Regression Testing, with a Bioinformatics Case Study.** Submitted for publication.
> Preprint DOI: *to be added once minted*.

A companion technical note derived from part of this study's positive-control panel is
published separately:

> Mahjoubi H. **False Novelty from Taxid Retagging in Historical NCBI Taxonomy
> Reconstruction.** Preprint: Zenodo, https://doi.org/10.5281/zenodo.21861803. Code:
> https://github.com/habibmahjoubi/taxid-retagging-artifact (not duplicated here).

## Summary

This repository supports a case study testing whether the magnitude of a knowledge-base
change (a real viral-taxonomy revision) predicts the resulting change in a downstream
classification pipeline's output -- the empirical basis for the paper's Knowledge-Drift-
Calibrated Revalidation (KDCR) meta-model. Two spiked viral targets (RSV, REO) were
reclassified with Kraken2 against thirteen historical reference databases, one per real
ICTV Master Species List transition (2014-2025), alongside a genuinely novel taxon used as
a positive control. The main finding: once taxonomic rank and reference-library composition
are held constant, whole-taxonomy drift magnitude (delta_K) does not predict classification
change for either established target, while the positive control shows the pipeline can and
does detect a genuine composition change -- i.e. a structural diff's *size* does not scope
its *impact*; content-level analysis is needed.

## Repository contents

```
scripts/    All analysis, database-construction, classification, and statistics scripts
            (see table below); shell scripts assume ~/kdcr as the base working directory.
results/    Summary CSV/TXT tables consumed by the manuscript's tables and figures.
positive_control/             Full working set for the Canine distemper virus positive
                               control (simulated reads, Kraken2 reports/outputs).
positive_control_conf/        Confidence-threshold sensitivity variant of the above.
positive_control_robustness/  Sequencing-error-rate sensitivity variant of the above.
```

Excluded from this repository (regenerable from the public sources below via the scripts
here, but too large to version -- tens of GB): the raw downloaded sequencing reads, the
built Kraken2 databases themselves, the downloaded NCBI taxdump snapshots, and the raw
per-read Kraken2 classification outputs for the full 13-database x 9-run panel (their
extracted summary statistics are the CSVs kept under `results/`).

## Script index

| Script | Purpose |
|---|---|
| **Reference-database construction** | |
| `tag_viral_fasta_with_taxid.py` | Tag the bulk RefSeq viral FASTA with `\|kraken:taxid\|N` headers via batched NCBI E-utilities `esummary` calls. |
| `filter_fasta_by_taxid_existence.py` | Filter a tagged FASTA to taxids valid in a given historical `nodes.dmp` snapshot -- operationalizes "knowledge drift" as library composition. |
| `build_historical_dbs.sh` | Build one Kraken2 database per ICTV MSL version (identical sequence set, different historical taxonomy tree). |
| `build_fixed_library.py` / `build_fixed_dbs.sh` | Build the confound-controlled variant: a single library fixed to taxids valid in *every* snapshot, isolating structural/topology change from composition change. |
| **Data acquisition** | |
| `download_runs.sh` | Subsample (ENA byte-range partial fetch) the 9 case-study runs from BioProject PRJNA1026487. |
| `download_runs_large.sh` / `download_one_run_large.sh` | 5x-deeper (10M read-pair) resample of the four informative runs for a statistical-power check. |
| **Classification** | |
| `classify_runs.sh` | Classify all 9 runs against the current reference database (baseline). |
| `classify_all_versions.sh` | Classify all 9 runs against each per-MSL-version database. |
| `classify_fixed_dbs.sh` | Classify against the confound-controlled fixed-library databases. |
| `classify_large_runs.sh` | Classify the 5x-deeper resamples against the fixed-library databases. |
| `classify_conf_sweep.sh` | Repeat classification at non-default `--confidence` thresholds. |
| **Knowledge-drift (delta_K) quantification** | |
| `compute_delta_K_extended.py` | Species-level delta_K across all twelve real MSL transitions (2014-2025). |
| `compute_delta_K_ncbi.py` | delta_K computed directly from consecutive NCBI taxdump snapshots (same source that drives classification). |
| `compute_delta_K_local_formal.py` / `compute_candidate1_rankmatched.py` | Formal delta_K_local candidates from the NCBI taxonomy graph (node-relative and descendant-drift). |
| **Core result summaries** | |
| `summarize_classification.py` | RSV/REO read counts per run, current database. |
| `summarize_all_versions.py` / `summarize_all_versions_extended.py` | Read counts across all historical databases; change vs. delta_K. |
| `recompute_rank_matched.py` | Recompute REO at a rank-matched target taxid (removes the genus-vs-species rank confound). |
| `summarize_fixed.py` | Same summary against the confound-controlled fixed-library databases. |
| **Positive control (genuinely novel taxon)** | |
| `find_newly_added_taxon.py` | Identify a taxon whose taxid is absent from an older snapshot but present in a newer one. |
| `test_positive_control.sh` | Simulate and classify reads from the positive-control genome across the bracketing MSL databases. |
| `m1_quantify_merge_fix.py` / `m1_artifact_verify.py` | Merge-aware audit distinguishing genuinely novel taxa from taxid-renumbering artifacts in the 56-taxon panel (full detail in the companion technical note above). |
| `m7_coverage_matched_56taxa.py` / `stratify_56taxa.py` | Coverage-matched resimulation and stratified (phage vs. non-phage) analysis of the 56-taxon panel. |
| **Structural-edit (reparenting) panel** | |
| `find_reparented_fixed.py` | Identify fixed-library taxa that underwent genuine reparenting on the same transition as the positive control. |
| `reparenting_panel_test.py` | Simulate and classify reads for that panel -- tests whether structural editing alone (no composition change) moves classification. |
| **Robustness and statistical checks** | |
| `m5_permutation_joint_null.py` / `m5b_circular_permutation.py` | Family-wise and temporally-aware permutation nulls for the delta_K correlation tests. |
| `perm_test_ncbi_viral.py` | Permutation test for the delta_K_ncbi_viral correlation. |
| `analyze_large_runs.py` / `analyze_conf05_sweep.py` | Analyze the 5x-depth and `--confidence 0.5` robustness variants. |
| `check_gained_35_36.py` | Diagnostic check on a specific transition's library-composition change. |
| `read_concordance.py` / `characterize_discordance.py` | Per-read concordance analysis across MSL versions and characterization of the largest discordance spikes (Microviridae/PhiX174). |

## Dependencies

- Python 3.10+ with `pandas`, `numpy`, `scipy`, `matplotlib`, `biopython` (see `requirements.txt`)
- [Kraken2](https://github.com/DerrickWood/kraken2) (tested with v2.17.1) and `krakentools`
- [sra-tools](https://github.com/ncbi/sra-tools) and NCBI `entrez-direct`, on `PATH`
- [wgsim](https://github.com/lh3/wgsim) (or an equivalent read simulator) on `PATH`
- Internet access to NCBI E-utilities and the ENA/SRA mirrors, for data acquisition scripts

Set up the Kraken2/bioinformatics tooling (conda/mamba):

```
conda create -n kdcr -c bioconda -c conda-forge kraken2 krakentools sra-tools entrez-direct
conda activate kdcr
pip install -r requirements.txt
```

## Data

This study uses exclusively public data, none of which is redistributed in this repository:

- ICTV Master Species Lists MSL29-MSL41: https://ictv.global/msl
- NCBI taxonomy archive snapshots: https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump_archive/
- NCBI RefSeq viral genome release: https://ftp.ncbi.nlm.nih.gov/refseq/release/viral/
- Sequencing data: BioProject PRJNA1026487, via NCBI SRA and its ENA mirror

## Reproducing the results

Shell scripts assume a working directory `~/kdcr` with subdirectories `data/`,
`kraken2_dbs/`, `kraken2_dbs_fixed/`, `results/`, `results_fixed/`, matching the paths
referenced in the Python scripts. Broad sequence, per manuscript Methods section order:

1. Download the taxonomy snapshots and RefSeq viral release; tag and filter as needed
   (`tag_viral_fasta_with_taxid.py`, `filter_fasta_by_taxid_existence.py`).
2. Build the per-version and fixed-library databases (`build_historical_dbs.sh`,
   `build_fixed_library.py`, `build_fixed_dbs.sh`).
3. Download and subsample the case-study reads (`download_runs.sh`).
4. Classify (`classify_*.sh`) and summarize (`summarize_*.py`).
5. Compute delta_K (`compute_delta_K_*.py`) and correlate against classification change.
6. Run the positive-control and reparenting-panel scripts for the two contrast experiments.
7. Run the robustness/statistical-null scripts for the checks reported in the Supplementary
   Results and Supplementary Statistical Notes.

## License

Code in this repository is released under the MIT License (see `LICENSE`). The underlying
NCBI taxonomy, RefSeq, and SRA data are public and are not redistributed here.

## Citation

If you use this code, please cite the manuscript above. A permanent DOI for this exact
repository snapshot is provided by Zenodo (link to be added once minted).
