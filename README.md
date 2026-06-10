# Reciprocal I2 relic retention

Public replication package for:

**Reciprocal relic retention in the Serbian DNA Project Y-chromosome register: evidence for a broader I2 patrilineal reservoir**

The repository rebuilds the reported analyses from the public source tables and regenerates the result tables and figures. It contains only data, R code, result tables, figures, and repository metadata.

## Repository Contents

Core runner and metadata:

- `run_all.R`
- `config_analysis.R`
- `r-requirements.txt`
- `CITATION.cff`
- `LICENSE`
- `DATA_SOURCE.md`

Public source inputs:

- `dnk_source_snapshot.csv`
- `dnk_source_table.csv`
- `branch_ages.csv`
- `geocoded_locations_cache.rds`
- `ph908_branch_map.tsv`
- `ph908_node_age_tiers.tsv`

Numbered analysis scripts:

- `analysis/scripts/01_preprocess_registry_data.R`
- `analysis/scripts/02_build_branch_reference_table.R`
- `analysis/scripts/03_build_class_specific_units.R`
- `analysis/scripts/04_build_mutual_i2_units.R`
- `analysis/scripts/05_confirm_mutual_i2_coretention.R`
- `analysis/scripts/06_confirm_ph908_r1a_polarity.R`
- `analysis/scripts/07_exact_region_polarity_check.R`
- `analysis/scripts/08_geography_only_reciprocal_i2_test.R`
- `analysis/scripts/09_dual_channel_reciprocal_i2_test.R`
- `analysis/scripts/10_coordinate_surname_slava_ecology.R`
- `analysis/scripts/11_cross_surname_location_i2_coretention.R`
- `analysis/scripts/12_own_family_reciprocal_reservoir_test.R`
- `analysis/scripts/13_own_family_frequency_conditioned_nulls.R`
- `analysis/scripts/14_summarise_evidence_summary.R`

Figure and sensitivity scripts:

- `scripts/01_make_wider_i2_figures.R`
- `scripts/02_make_ph908_removed_sensitivity.R`
- `scripts/09_make_registry_bias_sensitivity.R`

Public-release outputs:

- `result_tables/`
- `figures/`

Generated intermediate folders such as `analysis/outputs/`, `results/`, `tables/`, `facts/`, `figuredata/`, and `logs/` are reproducible or local-only and are excluded from the lean public archive.

## Rebuild

Install R dependencies:

```powershell
Rscript -e "install.packages(readLines('r-requirements.txt'))"
```

Run the full analysis pipeline from the repository root:

```powershell
Rscript .\run_all.R full
```

The default mode is also `full`:

```powershell
Rscript .\run_all.R
```

Full mode uses 5,000 permutations for the reported confirmation tests. It rebuilds preprocessing outputs, analysis result tables, figures, and `pipeline_manifest.csv`.

A quick validation mode is available:

```powershell
Rscript .\run_all.R quick
```

Quick mode is only for checking that the pipeline executes in a local environment. It does not reproduce the final reported p-values.

Figures-only mode rebuilds the figures and manifest from existing `result_tables/`:

```powershell
Rscript .\run_all.R figures
```

## Analysis Scope

The analysis asks whether the broader I2 lineage family shows reciprocal relic retention in the Serbian DNA Project Y-chromosome register. Each sufficiently represented Y-chromosome family retained for analysis is scored against its own relic-class field using the same rule. Same-surname and same-location echoes are excluded before scoring. The package includes social, coordinate, hard frequency-conditioned, deletion, PH908-removed, PH908/R1a polarity, and registry-bias sensitivity tests.

The primary reported result is the surname-region-terminal own-family social difference-in-differences statistic. The coordinate-cloud test is an independent spatial companion. The family ranking, PH908/R1a polarity bridge, PH908-removed dominance sensitivity, and registry-bias tipping-point analysis are robustness and interpretation layers.

## Main Outputs

- `result_tables/`
- `figures/`

## Citation

Until the Zenodo DOI is minted, cite the repository as:

Stanisic, N. *Reciprocal I2 Relic Retention Public Replication Package* (v1.0.0). GitHub repository: `https://github.com/nstanisic1/wider-i2-reciprocal-reservoir`.

## License

Code is released under the MIT License. Source data provenance is described in `DATA_SOURCE.md`.
