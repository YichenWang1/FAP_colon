---
title: "Mutational signature hypothesis testing"
author: "Yichen Wang"
output:
  html_document:
    keep_md: yes
---



## Setup


``` r
suppressPackageStartupMessages({
  library(nlme)
})
```


``` r
ROOT    <- "/lustre/scratch125/casm/teams/team267/users/yw2/colon/github_repo"

OUT_DIR <- file.path(ROOT, "Fig3_Signature_analysis")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

FAP_CSV  <- file.path(ROOT, "data/fap_database_withexposure.csv")
CTRL_CSV <- file.path(ROOT, "data/ctrl_database_withexposure.csv")

# Lesion stages in order of progression; each stage is compared with normal
# crypts in a model that also contains every earlier stage.
LADDER    <- c("normal", "aca", "microadenoma", "diminutive_polyp", "large_polyp")
# Signatures tested along the ladder
LADDER_SIGS <- c("SBS18", "SBS1", "SBS5", "ID1", "ID2", "ID5")
# Signatures whose per-year rate is compared between FAP and wild-type normals
BASELINE_SIGS <- c("SBS1", "SBS5", "SBS18", "ID1", "ID2", "ID5")
# Ladder random effects: subgroup_id separates independent APC-mutant lineages
# within a lesion (matches burden_hypothesis_testing and Fig3bc_disease_effect).
LADDER_RANDOM <- ~1 | patient/block/subgroup_id
# Normal-crypt-only models have no lesion structure below the biopsy block.
NORMAL_RANDOM <- ~1 | patient/block
```

## Build the modelling dataset


``` r
database <- read.csv(FAP_CSV, header = TRUE, stringsAsFactors = FALSE,
                       na.strings = c("", "NA"))
# monocryptal and oligocryptal adenomas are analysed together as microadenomas
database$type[database$type %in% c("monocryptal_adenoma", "oligocryptal_adenoma")] <- "microadenoma"
database$type <- relevel(as.factor(database$type), ref = "normal")

# germline + somatic APC hits: 1 = APC mut/wt, >=2 = second hit present
database$n_apc_hits <- as.numeric(3 - (is.na(database$germlinehit) +
                                        is.na(database$somatic1) +
                                        is.na(database$somatic2)))
database <- database[database$coverage > 10, ]
database$germlinehit <- as.factor("APC")

# Normal crypts carry the literal group_id/subgroup_id 'normal', which would
# make them one clonal unit in the nested random effects. Give each its own id.
is_normal <- database$type == "normal"
database$group_id[is_normal]    <- database$sample[is_normal]
database$subgroup_id[is_normal] <- database$sample[is_normal]

# Non-adenomatous lesion crypts without a second APC hit are bystanders, not
# lesion crypts; they get their own level and so drop out of the ladder.
database$type <- factor(database$type,
                          levels = c(levels(database$type), "bystander"))
database$type[database$n_apc_hits < 2 &
                database$type %in% c("rc", "bc", "aca", "trapped_crypt")] <- "bystander"

# Adenomatous crypts without a second APC hit are dropped (none at present).
adenoma_no_second_hit <- database$n_apc_hits < 2 &
  database$type %in% c("microadenoma", "diminutive_polyp", "large_polyp")
cat("dropped, adenomatous lesion crypts with no somatic APC mutation:",
    sum(adenoma_no_second_hit), "\n")
```

```
## dropped, adenomatous lesion crypts with no somatic APC mutation: 4
```

``` r
if (any(adenoma_no_second_hit)) print(table(droplevels(database$type[adenoma_no_second_hit])))
```

```
## 
## diminutive_polyp 
##                4
```

``` r
database <- database[!adenoma_no_second_hit, ]
```

## Whether the rate of each signature per year differs between FAP and wild-type


``` r
# Make a combined dataset of FAP and wild-type normals, for the age-slope comparison.
# germlinehit recoded as the cohort: "APC" (FAP/AFAP) or "WT" (reference).
fap_and_ctrl <- local({
  database      <- read.csv(FAP_CSV,  header = TRUE, stringsAsFactors = FALSE, na.strings = c("", "NA"))
  ctrl_database <- read.csv(CTRL_CSV, header = TRUE, stringsAsFactors = FALSE, na.strings = c("", "NA"))
  database      <- database[database$coverage > 10, ]
  ctrl_database <- ctrl_database[ctrl_database$coverage > 10, ]
  shared_cols <- intersect(colnames(database), colnames(ctrl_database))
  both <- rbind(database[, shared_cols], ctrl_database[, shared_cols])
  both$germlinehit[is.na(both$germlinehit)]  <- "WT"
  both$germlinehit[both$germlinehit != "WT"] <- "APC"
  both$germlinehit <- relevel(as.factor(both$germlinehit), ref = "WT")
  both
})

# Likelihood-ratio test per signature: separate age slopes for FAP and WT
# (sig ~ age:germlinehit) against one shared slope (sig ~ age), normal crypts only.
lrt_by_signature <- function(crypts, label) {
  normals <- crypts[crypts$type == "normal", ]
  cat("\n==", label, "--", sum(normals$germlinehit == "APC"), "FAP /",
      sum(normals$germlinehit == "WT"), "control normal crypts\n")
  out <- do.call(rbind, lapply(BASELINE_SIGS, function(sig) {
    fit_data <- normals[!is.na(normals[[sig]]), ]
    separate_slopes <- lme(as.formula(paste(sig, "~ age:germlinehit")),
                           random = NORMAL_RANDOM, data = fit_data, method = "ML")
    shared_slope    <- lme(as.formula(paste(sig, "~ age")),
                           random = NORMAL_RANDOM, data = fit_data, method = "ML")
    lrt  <- anova(separate_slopes, shared_slope)
    coef <- fixef(separate_slopes)
    data.frame(signature = sig,
               rate_WT   = round(coef[["age:germlinehitWT"]],  3),
               rate_FAP  = round(coef[["age:germlinehitAPC"]], 3),
               L.Ratio   = round(lrt[2, "L.Ratio"], 3),
               p         = signif(lrt[2, "p-value"], 3))
  }))
  out$p_adj <- signif(p.adjust(out$p, method = "BH"), 3)
  print(out, row.names = FALSE)
  invisible(cbind(comparison = label, out))
}

rate_results <- lrt_by_signature(fap_and_ctrl, "all FAP normal")
```

```
## 
## == all FAP normal -- 192 FAP / 400 control normal crypts
##  signature rate_WT rate_FAP L.Ratio        p   p_adj
##       SBS1  13.255   17.370   2.746 0.097500 0.14600
##       SBS5  18.007   23.509   5.556 0.018400 0.05520
##      SBS18   6.323    7.612   0.592 0.441000 0.44100
##        ID1   0.398    0.637   3.320 0.068400 0.13700
##        ID2   0.152    0.218   0.877 0.349000 0.41900
##        ID5   0.669    1.312  13.315 0.000263 0.00158
```

``` r
# normal crypts of both cohorts, used by the sequencing-depth check below
baseline_normals <- fap_and_ctrl[fap_and_ctrl$type == "normal", ]
```

### Is the higher ID5/borderline SBS5 rate in FAP normal crypts real? Controlling sequencing depth

The FAP crypts were sequenced deeper than the wild-type crypts (Lee-Six et al. 2019).
Flat signatures (SBS5, ID5) are where extra low-VAF calls tend to be assigned, so the
rate comparison is repeated (i) with coverage as a covariate and (ii) restricted to
crypts at >=20x.


``` r
cat("coverage quantiles:\n")
```

```
## coverage quantiles:
```

``` r
print(do.call(rbind, tapply(baseline_normals$coverage, baseline_normals$germlinehit,
                            quantile, c(0, .25, .5, .75, 1))))
```

```
##        0%   25%   50%     75%  100%
## WT  10.17 14.53 16.32 18.0725 34.55
## APC 11.45 26.25 31.10 35.5750 50.45
```

``` r
# One LRT: `separate_formula` (FAP/WT slopes) against `shared_formula`.
fit_rate_lrt <- function(crypts, separate_formula, shared_formula) {
  separate_slopes <- lme(separate_formula, random = NORMAL_RANDOM, data = crypts, method = "ML")
  shared_slope    <- lme(shared_formula,   random = NORMAL_RANDOM, data = crypts, method = "ML")
  coef <- fixef(separate_slopes)
  data.frame(rate_WT  = round(coef[["age:germlinehitWT"]], 2),
             rate_FAP = round(coef[["age:germlinehitAPC"]], 2),
             p = signif(anova(separate_slopes, shared_slope)[2, "p-value"], 3))
}

# Both depth controls for total SBS, total indel and each baseline signature.
depth_controlled_lrts <- function(normals, label) {
  deep <- normals[normals$coverage >= 20, ]
  cat(sprintf("\n%s -- crypts >=20x: FAP %d crypts / %d individuals; WT %d / %d\n", label,
              sum(deep$germlinehit == "APC"), length(unique(deep$patient[deep$germlinehit == "APC"])),
              sum(deep$germlinehit == "WT"),  length(unique(deep$patient[deep$germlinehit == "WT"]))))
  do.call(rbind, lapply(c("sbstotal_corr", "indeltotal_corr", BASELINE_SIGS), function(sig) {
    with_response <- function(rhs) as.formula(paste(sig, rhs))
    cbind(set = label, signature = sig, rbind(
      cbind(model = "+ coverage",
            fit_rate_lrt(normals, with_response("~ age:germlinehit + coverage"),
                                  with_response("~ age + coverage"))),
      cbind(model = "crypts >=20x",
            fit_rate_lrt(deep, with_response("~ age:germlinehit"), with_response("~ age")))))
  }))
}
res_depth <- depth_controlled_lrts(baseline_normals, "all FAP normal")
```

```
## 
## all FAP normal -- crypts >=20x: FAP 176 crypts / 15 individuals; WT 41 / 9
```

``` r
print(res_depth, row.names = FALSE)
```

```
##             set       signature        model rate_WT rate_FAP     p
##  all FAP normal   sbstotal_corr   + coverage   35.62    34.64 0.874
##  all FAP normal   sbstotal_corr crypts >=20x   46.14    45.93 0.976
##  all FAP normal indeltotal_corr   + coverage    1.13     0.79 0.462
##  all FAP normal indeltotal_corr crypts >=20x    2.44     2.59 0.831
##  all FAP normal            SBS1   + coverage   13.65    15.96 0.357
##  all FAP normal            SBS1 crypts >=20x   17.17    17.56 0.896
##  all FAP normal            SBS5   + coverage   18.52    21.38 0.211
##  all FAP normal            SBS5 crypts >=20x   20.04    20.90 0.754
##  all FAP normal           SBS18   + coverage    6.39     7.38 0.565
##  all FAP normal           SBS18 crypts >=20x    6.96     7.29 0.746
##  all FAP normal             ID1   + coverage    0.46     0.41 0.703
##  all FAP normal             ID1 crypts >=20x    0.95     0.96 0.960
##  all FAP normal             ID2   + coverage    0.18     0.14 0.613
##  all FAP normal             ID2 crypts >=20x    0.47     0.43 0.734
##  all FAP normal             ID5   + coverage    0.78     0.90 0.399
##  all FAP normal             ID5 crypts >=20x    1.07     1.42 0.174
```

## The ladder

`normal-> aca -> microadenoma -> diminutive_polyp -> large_polyp`, each
category fitted against the categories at or below it. All model output is echoed here, in
the knitted document.


``` r
invisible(capture.output({
  stage_results <- list()
  for (k in 2:length(LADDER)) {
    stage <- LADDER[k]
    stages_in_model <- LADDER[1:k]
    cat("\n# Whether", stage, "leads to increased burden\n\n")
    stage_rows <- list()
    for (sig in LADDER_SIGS) {
      fit_data <- database[database$type %in% stages_in_model, ]
      fit_data$type <- factor(fit_data$type, levels = stages_in_model)
      # REML: these fits are reported (slope, SE, CI, Wald p), never LR-tested
      fit <- lme(as.formula(paste0(sig, " ~ age + type")),
                 random = LADDER_RANDOM, data = fit_data, method = "REML")
      cat("--", sig, ",", stage, "--\n")
      print(summary(fit))
      print(intervals(fit, which = "fixed"))
      # the stage's own coefficient: its excess burden over normal at a given age
      coef_name <- paste0("type", stage)
      coef_row <- summary(fit)$tTable[coef_name, ]
      ci       <- intervals(fit, which = "fixed")$fixed[coef_name, ]
      stage_rows[[sig]] <- data.frame(
        response = sig, estimate = coef_row[["Value"]], std_error = coef_row[["Std.Error"]],
        ci_lower = ci[["lower"]], ci_upper = ci[["upper"]], p = coef_row[["p-value"]],
        p_adj = NA_real_)
    }
    stage_table <- do.call(rbind, stage_rows)
    rownames(stage_table) <- NULL
    cat("\n## raw P table,", stage, "(BH adjustment is applied across the whole ladder below)\n")
    print(stage_table[, setdiff(names(stage_table), "p_adj")])
    stage_results[[stage]] <- stage_table
  }
}))
```

## BH-adjusted results


``` r
combined <- do.call(rbind, lapply(names(stage_results),
                                  function(stage) cbind(category = stage, stage_results[[stage]])))
# one family: all signature x category tests
combined$p_adj <- p.adjust(combined$p, method = "BH")
cat("BH-adjusted across all", nrow(combined), "signature x category tests\n")
```

```
## BH-adjusted across all 24 signature x category tests
```

``` r
results <- data.frame(
  comparison = "ladder", group = combined$category, response = combined$response,
  estimate = combined$estimate, std_error = combined$std_error,
  ci_lower = combined$ci_lower, ci_upper = combined$ci_upper,
  L.Ratio = NA_real_, p = combined$p, p_adj = combined$p_adj,
  stringsAsFactors = FALSE)

out_csv <- file.path(OUT_DIR, "signature_results.csv")
write.csv(results, out_csv, row.names = FALSE)

knitr::kable(combined, digits = 4, row.names = FALSE)
```



|category         |response |  estimate| std_error|  ci_lower|  ci_upper|      p|  p_adj|
|:----------------|:--------|---------:|---------:|---------:|---------:|------:|------:|
|aca              |SBS18    |   95.8351|   20.3443|   55.6735|  135.9967| 0.0000| 0.0000|
|aca              |SBS1     |  138.7722|   34.1018|   71.4518|  206.0926| 0.0001| 0.0001|
|aca              |SBS5     |   34.5370|   25.7012|  -16.1998|   85.2738| 0.1808| 0.1887|
|aca              |ID1      |   21.6353|    2.6104|   16.4821|   26.7886| 0.0000| 0.0000|
|aca              |ID2      |    8.4274|    1.3468|    5.7686|   11.0862| 0.0000| 0.0000|
|aca              |ID5      |    1.5427|    3.5353|   -5.4365|    8.5218| 0.6631| 0.6631|
|microadenoma     |SBS18    |  333.2792|   28.6105|  276.8176|  389.7409| 0.0000| 0.0000|
|microadenoma     |SBS1     |  273.2221|   36.5387|  201.1145|  345.3296| 0.0000| 0.0000|
|microadenoma     |SBS5     |  127.7571|   25.5609|   77.3138|  178.2004| 0.0000| 0.0000|
|microadenoma     |ID1      |   34.4323|    3.0433|   28.4265|   40.4380| 0.0000| 0.0000|
|microadenoma     |ID2      |   20.6623|    1.7509|   17.2071|   24.1176| 0.0000| 0.0000|
|microadenoma     |ID5      |   12.2303|    3.3080|    5.7021|   18.7585| 0.0003| 0.0003|
|diminutive_polyp |SBS18    |  530.9431|   34.7802|  462.3472|  599.5390| 0.0000| 0.0000|
|diminutive_polyp |SBS1     |  504.9193|   35.5789|  434.7481|  575.0905| 0.0000| 0.0000|
|diminutive_polyp |SBS5     |  325.8080|   41.2343|  244.4828|  407.1331| 0.0000| 0.0000|
|diminutive_polyp |ID1      |   74.4686|    5.9086|   62.8152|   86.1220| 0.0000| 0.0000|
|diminutive_polyp |ID2      |   39.4476|    5.0270|   29.5329|   49.3623| 0.0000| 0.0000|
|diminutive_polyp |ID5      |   25.6064|    2.6385|   20.4025|   30.8103| 0.0000| 0.0000|
|large_polyp      |SBS18    | 1751.8579|   96.3587| 1561.8128| 1941.9030| 0.0000| 0.0000|
|large_polyp      |SBS1     | 1304.7840|   91.5520| 1124.2190| 1485.3490| 0.0000| 0.0000|
|large_polyp      |SBS5     | 1474.6475|  137.3212| 1203.8134| 1745.4816| 0.0000| 0.0000|
|large_polyp      |ID1      |  395.4544|   20.1891|  355.6361|  435.2728| 0.0000| 0.0000|
|large_polyp      |ID2      |  258.5005|   14.8860|  229.1413|  287.8598| 0.0000| 0.0000|
|large_polyp      |ID5      |   76.1637|    6.9801|   62.3971|   89.9304| 0.0000| 0.0000|


``` r
significant <- combined[combined$p_adj < 0.05, ]
knitr::kable(significant[order(significant$category, significant$p_adj), ],
             digits = 4, row.names = FALSE,
             caption = "Signature/category combinations reaching P-adjust < 0.05")
```



Table: Signature/category combinations reaching P-adjust < 0.05

|category         |response |  estimate| std_error|  ci_lower|  ci_upper|     p| p_adj|
|:----------------|:--------|---------:|---------:|---------:|---------:|-----:|-----:|
|aca              |ID1      |   21.6353|    2.6104|   16.4821|   26.7886| 0e+00| 0e+00|
|aca              |ID2      |    8.4274|    1.3468|    5.7686|   11.0862| 0e+00| 0e+00|
|aca              |SBS18    |   95.8351|   20.3443|   55.6735|  135.9967| 0e+00| 0e+00|
|aca              |SBS1     |  138.7722|   34.1018|   71.4518|  206.0926| 1e-04| 1e-04|
|diminutive_polyp |SBS18    |  530.9431|   34.7802|  462.3472|  599.5390| 0e+00| 0e+00|
|diminutive_polyp |SBS1     |  504.9193|   35.5789|  434.7481|  575.0905| 0e+00| 0e+00|
|diminutive_polyp |ID1      |   74.4686|    5.9086|   62.8152|   86.1220| 0e+00| 0e+00|
|diminutive_polyp |ID5      |   25.6064|    2.6385|   20.4025|   30.8103| 0e+00| 0e+00|
|diminutive_polyp |SBS5     |  325.8080|   41.2343|  244.4828|  407.1331| 0e+00| 0e+00|
|diminutive_polyp |ID2      |   39.4476|    5.0270|   29.5329|   49.3623| 0e+00| 0e+00|
|large_polyp      |ID1      |  395.4544|   20.1891|  355.6361|  435.2728| 0e+00| 0e+00|
|large_polyp      |SBS18    | 1751.8579|   96.3587| 1561.8128| 1941.9030| 0e+00| 0e+00|
|large_polyp      |ID2      |  258.5005|   14.8860|  229.1413|  287.8598| 0e+00| 0e+00|
|large_polyp      |SBS1     | 1304.7840|   91.5520| 1124.2190| 1485.3490| 0e+00| 0e+00|
|large_polyp      |ID5      |   76.1637|    6.9801|   62.3971|   89.9304| 0e+00| 0e+00|
|large_polyp      |SBS5     | 1474.6475|  137.3212| 1203.8134| 1745.4816| 0e+00| 0e+00|
|microadenoma     |ID2      |   20.6623|    1.7509|   17.2071|   24.1176| 0e+00| 0e+00|
|microadenoma     |SBS18    |  333.2792|   28.6105|  276.8176|  389.7409| 0e+00| 0e+00|
|microadenoma     |ID1      |   34.4323|    3.0433|   28.4265|   40.4380| 0e+00| 0e+00|
|microadenoma     |SBS1     |  273.2221|   36.5387|  201.1145|  345.3296| 0e+00| 0e+00|
|microadenoma     |SBS5     |  127.7571|   25.5609|   77.3138|  178.2004| 0e+00| 0e+00|
|microadenoma     |ID5      |   12.2303|    3.3080|    5.7021|   18.7585| 3e-04| 3e-04|
