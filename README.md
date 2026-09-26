# kdcr-content-vs-structure
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.21878420.svg)](https://doi.org/10.5281/zenodo.21878420)

Scripts, result tables and small input files for the article:

> Mahjoubi H. **Viral taxonomy updates rarely change Kraken2 detection of adventitious viruses.**
> Manuscript submitted for publication.

Companion article on taxonomic identifier artefacts (code in its own repository):

> Mahjoubi H. **When current taxids make historical taxa look new.** J. Bioinform. Comput. Biol.
> 24, 2671004 (2026). https://doi.org/10.1142/S0219720026710046 ·
> https://github.com/habibmahjoubi/taxid-retagging-artifact

## Study in brief

Twelve annual updates of the NCBI viral taxonomy (ICTV MSL29 to MSL41, 2014–2025) were replayed with
Kraken2 on public HTS runs of CHO-K1 harvests spiked with respiratory syncytial virus or mammalian
orthoreovirus (BioProject PRJNA1026487), on reads simulated from fourteen adventitious agents and from
the spiked orthoreovirus strain, with three database designs and 17 dated prebuilt Kraken2 viral
databases (2020–2026):

- **Design 1**, taxonomy-filtered library: a sequence is kept if its current taxid exists in the snapshot;
- **Design 2**, fixed, incomplete library: only taxids valid in all 13 snapshots;
- **Design 3**, fixed, complete library: all sequences, each labelled with the taxid valid at the date.

## Repository layout

```
scripts/                      all scripts (index below); shell scripts assume the working directory ~/kdcr
results/                      result tables (index below)
sequences/                    genomes simulated from outside RefSeq, 7SL RNA reference, reads aligned to NCBI nt
positive_control/             canine distemper virus reads classified with the genome absent and present
positive_control_conf/        same, at confidence thresholds 0.1 and 0.5
positive_control_robustness/  same, at substitution error rates of 2%, 5% and 8%
position_bias_check/          read-position bias check of the subsampled runs
reparenting_panel/            manifest of the first reparenting panel
```

Not versioned (tens of GB, rebuilt from public data with the scripts): sequencing reads, Kraken2
databases, NCBI taxonomy snapshots and per-read Kraken2 outputs.

## Scripts

| Step | Scripts |
|---|---|
| Reference databases | `tag_viral_fasta_with_taxid.py`, `filter_fasta_by_taxid_existence.py`, `build_historical_dbs.sh` (Design 1), `build_fixed_library.py`, `build_fixed_dbs.sh` (Design 2), `design3_labels.py` (Design 3) |
| Reads | `download_runs.sh`, `download_runs_large.sh`, `download_one_run_large.sh` (ENA subsamples); `agent_panel.py` (fourteen-agent panel); `reparent_all_updates.py`, `find_reparented_fixed.py`, `reparenting_panel_test.py` (reparenting); `find_newly_added_taxon.py`, `test_positive_control.sh`, `coverage_matched_56taxa.py` (missing-genome panel) |
| Classification | `classify_runs.sh`, `classify_all_versions.sh`, `classify_fixed_dbs.sh`, `classify_large_runs.sh`, `classify_conf_sweep.sh`, `classify_SRR26352204.sh`; `pipeline_designs_panel.sh` (Design 3, panel, reparenting, factorial builds) and `pipeline_exposure_host.sh` (exposure, host reads, prebuilt databases) |
| Target counts and exposure | `summarize_classification.py`, `summarize_all_versions.py`, `summarize_all_versions_extended.py`, `summarize_fixed.py`, `recompute_rank_matched.py`, `summarize_designs_panel.py`, `exposure_host_prebuilt.py`, `reparent_summary.py`, `reparent_lineage_exposure.py` |
| Read-level reassignment | `read_concordance.py`, `characterize_discordance.py`, `characterise_reassignment.py`, `reassignment_denominators.py`, `reassignment_merge_corrected.py` |
| Alignment of low counts and host reads | `alignment_and_simulation_checks.py`, `align_target_node_reads.sh`, `spurious_mrv_longest_match.py`, `reads_lost_at_confidence01.py`, `reference_record_contamination.sh` |
| Magnitude of taxonomic change | `compute_delta_K_extended.py` (species level), `compute_delta_K_ncbi.py` (from NCBI snapshots) |
| Missing-genome panel | `merge_aware_library_sizes.py`, `merge_artifact_verification.py`, `stratify_56taxa.py` |
| Statistical checks | `permutation_joint_null.py`, `permutation_circular.py`, `perm_test_ncbi_viral.py`, `compute_delta_K_local_formal.py`, `compute_candidate1_rankmatched.py`, `analyze_large_runs.py`, `analyze_conf05_sweep.py`, `check_gained_35_36.py` |

## Results

| Analysis | Files in `results/` |
|---|---|
| Magnitude of taxonomic change per update | `delta_K_species_extended_MSL30_41.csv`, `delta_K_ncbi_transitions.csv`, `delta_K_ncbi_vs_ictv_merged.csv` |
| Reads at the RSV and MRV target nodes in the spiked runs, all databases and confidence thresholds | `spiked_counts.tsv`, `confidence_counts.tsv`, `classification_*.csv`, `table3_fixed_library_raw_counts.csv`, `rsv_reo_summary_current_db.csv`, `current_db_conf/`, `cross_negatives.tsv`, `SRR26352204_all_designs.tsv` |
| Alignment of every target-node read (viral library, NCBI nt, Chinese hamster genome) | `target_reads_blast.tsv`, `target_reads_blast_nt.txt`, `reads_at_target_nodes_all_databases.tsv`, `reads_at_target_nodes_alignment.tsv` |
| Reads simulated from the targets and from the spiked MRV type 1 Lang strain | `simulated_targets_counts.tsv`, `simulated_targets_read_changes.tsv`, `mrv_type1_lang_simulated.tsv` |
| Panel of fourteen adventitious agents, exposure bound | `panel_manifest.tsv`, `panel_by_version.tsv`, `panel_changes.tsv`, `official_panel_by_version.tsv`, `official_panel_changes.tsv`, `panel_exposure_by_version.tsv`, `panel_exposure_summary.tsv`, `design3_panel_labels.tsv`, `static_sequences_under_target.tsv`, `vesivirus2117_*.tsv` |
| Reparenting at fixed content, all updates | `reparent_manifest.tsv`, `reparent_skipped_short.tsv`, `reparent_all_updates.tsv`, `reparent_accessions_changing.tsv`, `reparent_lineage_exposure.tsv` |
| Genome missing from the reference (canine distemper virus and 56-sequence panel) | `positive_control_56taxa.csv`, `positive_control_56taxa_covmatched.csv`, `positive_control_56taxa_stratified.csv`, `merge_artifact_verification.csv`, `p56_misassigned_destination_changes.tsv` |
| Read-level reassignment between versions (Design 2, prebuilt databases) | `reassignment_denominators.tsv`, `reassignment_without_phix.tsv`, `reassignment_alternative_pairing.tsv`, `reassignment_merge_corrected.tsv`, `reassigned_reads_origin.tsv`, `reassigned_top_taxid_pairs.tsv`, `phix_in_fixed_libraries.tsv`, `characterise_jumps.tsv`, `official_read_reassignment.tsv`, `classified_reads_per_version_9runs.tsv` |
| Prebuilt databases and non-viral sequence in reference records | `prebuilt_classified_totals.tsv`, `prebuilt_classified_jumps_top_taxa.tsv`, `prebuilt_confidence0.1_target_reads.tsv`, `rsv1e6_reads_conf0_vs_conf0.1_prebuilt.tsv` |
| Spurious assignments and database builds (factorial comparison, rebuild determinism) | `factorial_databases.tsv`, `factorial_mrv_clade_reads.tsv`, `factorial_overlap.tsv`, `build_determinism_14_runs.tsv`, `mrv_current_vs_msl41.tsv`, `spurious_mrv_longest_match.tsv` |
| Design 3 labelling per version | `design3_labelling_summary.tsv` |
| Statistical checks (correlations with the magnitude of change, permutation nulls) | `delta_K_local_formal.csv`, `classification_change_vs_delta_K*.csv`, `confidence_sensitivity/` |

## Dependencies

- Kraken2 v2.17.1, wgsim, BLAST+ (megablast), minimap2, sra-tools, NCBI entrez-direct
- Python 3.10+ with the packages in `requirements.txt`

```
conda create -n kdcr -c bioconda -c conda-forge kraken2 krakentools sra-tools entrez-direct blast minimap2 wgsim
conda activate kdcr
pip install -r requirements.txt
```

## Data

All data are public and are not redistributed here: ICTV Master Species Lists (https://ictv.global/msl);
NCBI Taxonomy snapshots (https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump_archive/); NCBI RefSeq viral
release (https://ftp.ncbi.nlm.nih.gov/refseq/release/viral/); sequencing reads, BioProject PRJNA1026487
(NCBI SRA / ENA); prebuilt Kraken2 viral indexes (https://benlangmead.github.io/aws-indexes/k2);
Chinese hamster genome assemblies GCF_003668045.3 and GCF_000223135.1.

## Reproducing the results

1. Build the reference libraries and databases (Reference databases row above).
2. Download and subsample the runs; simulate the panel, target and reparenting reads.
3. Run the classification scripts and the two pipelines (`pipeline_designs_panel.sh`, then
   `pipeline_exposure_host.sh`).
4. Run the summary and alignment scripts; each writes its tables to `results/`.

## License and citation

Code under the MIT License (`LICENSE`). Please cite the article above; the Zenodo concept DOI
[10.5281/zenodo.21878420](https://doi.org/10.5281/zenodo.21878420) resolves to the latest archived
version of this repository.
