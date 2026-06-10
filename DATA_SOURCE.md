# Data Source

The registry source tables in this repository were assembled from publicly available entries in the Srpski DNK projekat (Serbian DNA Project) public Y-DNA register.

The included source files are:

- `dnk_source_snapshot.csv`
- `dnk_source_table.csv`
- `branch_ages.csv`
- `geocoded_locations_cache.rds`
- `ph908_branch_map.tsv`
- `ph908_node_age_tiers.tsv`

`dnk_source_snapshot.csv` preserves the source snapshot used for provenance. `dnk_source_table.csv` is the cleaned source table consumed by preprocessing. The remaining files provide branch-age, geocoding, and PH908 branch-annotation inputs required to reproduce the analysis.

The repository also generates downstream analysis tables from these inputs through `run_all.R`. Generated tables should not be treated as independent source data.

Code in this repository is released under the MIT License. The source registry records remain attributed to the public source from which they were assembled.
