# Disinfectant-tolerance DNABERT2 feasibility summary

| model | balanced_accuracy_mean | balanced_accuracy_std | roc_auc_mean | roc_auc_std |
| --- | --- | --- | --- | --- |
| Known-QAC determinant BLAST rule | 0.967 |  | 0.967 |  |
| Best sensitivity strategy: 64 windows, mean, random_forest | 0.575 | 0.073 | 0.659 | 0.061 |
| DNABERT2 16-window mean + logistic regression | 0.616 | 0.102 | 0.638 | 0.111 |
| Capped k-mer TF-IDF | 0.46 | 0.051 | 0.521 | 0.075 |

## Interpretation

- Treat these results as a feasibility package, not a deployable model.
- The endpoint is low-level BC MIC tolerance, not use-level sanitizer survival.
- External validation remains the next requirement for a stronger manuscript claim.
