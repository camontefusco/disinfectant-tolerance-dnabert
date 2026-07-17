# Data Manifests

Create `isolate_manifest.csv` from the published supplementary tables after
checking the reuse terms for each dataset.

Required columns:

| Column | Meaning |
| --- | --- |
| `isolate_id` | Stable identifier for one bacterial isolate |
| `assembly_path` | Local FASTA assembly path |
| `label` | Binary endpoint: `1` tolerant, `0` sensitive |
| `group` | Clonal complex, PopPUNK cluster, or another lineage-aware group |

Recommended provenance columns:

| Column | Meaning |
| --- | --- |
| `mic_bc_mg_l` | Benzalkonium-chloride MIC in mg/L |
| `source` | Dataset or publication source |
| `ena_accession` | ENA run or assembly accession |

The initial Gmeiner et al. endpoint uses `mic_bc_mg_l >= 1.25` as the tolerant
class. Use-level sanitizer-survival outcomes from Harrand et al. are a separate
endpoint and should not be merged into the MIC label.

