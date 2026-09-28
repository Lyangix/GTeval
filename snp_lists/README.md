# SNP lists

Variant definitions used for the Height, BMI, and WHR gPRS scores. Each set has
an annotated CSV and a one-ID-per-line `.snplist` file with no header.

| Trait | GWS | Suggestive | Liberal |
| --- | ---: | ---: | ---: |
| Height | 2,313 | 16,020 | 30,937 |
| BMI | 435 | 4,788 | 13,980 |
| WHR | 172 | 1,991 | 7,485 |

GWS uses `P < 5e-8` and LD r-squared 0.1. Suggestive uses `P < 1e-5` and
LD r-squared 0.5. Liberal uses `P < 1e-3` and LD r-squared 0.5.
Sets were pruned separately.

## Files

Names follow `{trait}_{threshold}`, in lowercase. For example:

- [height_liberal.snplist](height_liberal.snplist): SNP IDs only.
- [height_liberal_gprs_snp_effect_allele_weights.csv](height_liberal_gprs_snp_effect_allele_weights.csv): variant annotations and weights.


## Use a list

From the repository root:

```r
variants <- read.csv("snp_lists/height_liberal_gprs_snp_effect_allele_weights.csv")
snp_ids <- readLines("snp_lists/height_liberal.snplist")
# For a local matrix whose column names match these identifiers:
# G_selected <- G[, match(snp_ids, colnames(G)), drop = FALSE]
```

The ID-only files also work as PLINK `--extract` lists when the variant IDs match.
