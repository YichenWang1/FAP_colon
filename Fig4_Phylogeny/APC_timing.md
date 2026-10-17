---
title: "APC timing"
author: "Yichen Wang"
output: 
  html_document:
    keep_md: yes
---





``` r
library(ggtree)
library(ape)
library(nlme)
library(dplyr)
```


``` r
# Paths via ROOT, as trunk_length.Rmd does: this file used one absolute path
# and two relative ones, so it only worked when knitted from its own directory.
ROOT <- "/lustre/scratch125/casm/teams/team267/users/yw2/colon"
ONSET_CSV <- file.path(ROOT, "github_repo/data/phylogenetic_tree/APC_onset_and_trunk_length.csv")

database = read.csv(file.path(ROOT, "metadata/fap_database_withexposure.csv"), header = T)
database$type = as.factor(database$type)
database$type = relevel(database$type, ref = "normal")

# add number of APC mutations
database$num_APC = 3 - (is.na(database$germlinehit) + is.na(database$somatic1) +
    is.na(database$somatic2))
database$num_APC = as.numeric(database$num_APC)

database = database[database$coverage > 10, ]
database$germlinehit = "APC"

database_exclude_bystander = database[(database$trapped_in_polyp != "" | database$num_APC >=
    2 | database$type == "normal"), ]
database_exclude_bystander$group_id[database_exclude_bystander$group_id == "normal"] = database_exclude_bystander$sample[database_exclude_bystander$group_id ==
    "normal"]
```


``` r
df = database_exclude_bystander[database_exclude_bystander$type %in% c("normal"),
    ]
lmm_sbs5 <- lme(SBS5 ~ age - 1, random = ~age - 1 | patient, data = df, method = "ML")
summary(lmm_sbs5)
```

```
## Linear mixed-effects model fit by maximum likelihood
##   Data: df 
##       AIC      BIC    logLik
##   2196.75 2206.523 -1095.375
## 
## Random effects:
##  Formula: ~age - 1 | patient
##              age Residual
## StdDev: 5.057993 63.16124
## 
## Fixed effects:  SBS5 ~ age - 1 
##        Value Std.Error DF  t-value p-value
## age 26.77034    1.3405 14 19.97042       0
## 
## Standardized Within-Group Residuals:
##         Min          Q1         Med          Q3         Max 
## -2.48732387 -0.62559124 -0.00267703  0.51304241  3.74327063 
## 
## Number of Observations: 192
## Number of Groups: 15
```


``` r
df_trunk_length = read.csv(ONSET_CSV)
# This chunk writes back to the file it reads, so re-running it on its own
# output used to fail: left_join then found SBS5_random_slope_patient already
# present, suffixed it .x/.y, and the mutate below could not see it. Drop the
# derived columns first so the step is idempotent and does not depend on
# trunk_length.Rmd having just rewritten the file.
df_trunk_length <- df_trunk_length[, setdiff(names(df_trunk_length), c("SBS5_random_slope_patient",
    "SBS5_full_slope", "age_onset_upperlim"))]
fixed_slope <- fixef(lmm_sbs5)["age"]
ranef_patient <- ranef(lmm_sbs5) %>%
    tibble::rownames_to_column("patient") %>%
    rename(SBS5_random_slope_patient = age)


df_trunk_length <- df_trunk_length %>%
    left_join(ranef_patient, by = "patient")
df_trunk_length <- df_trunk_length %>%
    mutate(SBS5_full_slope = fixed_slope + SBS5_random_slope_patient, age_onset_upperlim = SBS5_trunk/SBS5_full_slope) %>%
    mutate(age_onset_upperlim = if_else(age_onset_upperlim > age, age, age_onset_upperlim))

write.csv(df_trunk_length, ONSET_CSV, quote = F, row.names = F)
```


``` r
pvalues_lst = c()
df_trunk_length$trunk_driver_1 = T
df_trunk_length$trunk_driver_1[df_trunk_length$trunk_somatic_driver > 1] = F

lmm <- lme(SBS_total_branch ~ trunk_driver_1, random = ~1 | patient/lineage, data = df_trunk_length,
    method = "ML")
summary(lmm)  #SBS1
```

```
## Linear mixed-effects model fit by maximum likelihood
##   Data: df_trunk_length 
##        AIC      BIC    logLik
##   3268.288 3284.854 -1629.144
## 
## Random effects:
##  Formula: ~1 | patient
##         (Intercept)
## StdDev:   0.1181166
## 
##  Formula: ~1 | lineage %in% patient
##         (Intercept) Residual
## StdDev:    929.7517 568.7837
## 
## Fixed effects:  SBS_total_branch ~ trunk_driver_1 
##                        Value Std.Error  DF   t-value p-value
## (Intercept)        1778.2127  259.9290 159  6.841148   0.000
## trunk_driver_1TRUE -543.9723  316.7301  38 -1.717464   0.094
##  Correlation: 
##                    (Intr)
## trunk_driver_1TRUE -0.821
## 
## Standardized Within-Group Residuals:
##         Min          Q1         Med          Q3         Max 
## -2.78759461 -0.24291617 -0.02987252  0.17876996  6.51676032 
## 
## Number of Observations: 203
## Number of Groups: 
##              patient lineage %in% patient 
##                    5                   44
```

``` r
df_trunk = df_trunk_length[, c(1:5, 7:10, 17)]
df_trunk = unique(df_trunk)

lmm <- lme(SBS_total_trunk ~ trunk_driver_1, random = ~1 | patient, data = df_trunk,
    method = "ML")
summary(lmm)  #SBS5,SBStotal
```

```
## Linear mixed-effects model fit by maximum likelihood
##   Data: df_trunk 
##       AIC      BIC   logLik
##   794.172 801.3088 -393.086
## 
## Random effects:
##  Formula: ~1 | patient
##         (Intercept) Residual
## StdDev:    344.5982 1806.932
## 
## Fixed effects:  SBS_total_trunk ~ trunk_driver_1 
##                        Value Std.Error DF   t-value p-value
## (Intercept)         2557.156  525.3946 38  4.867116  0.0000
## trunk_driver_1TRUE -1523.601  603.7695 38 -2.523481  0.0159
##  Correlation: 
##                    (Intr)
## trunk_driver_1TRUE -0.784
## 
## Standardized Within-Group Residuals:
##        Min         Q1        Med         Q3        Max 
## -1.1977218 -0.5133457 -0.1261741  0.1749599  5.2436149 
## 
## Number of Observations: 44
## Number of Groups: 5
```

``` r
pvalues_lst = c(pvalues_lst, summary(lmm)$tTable[, "p-value"])
```

``` r
lmm <- lme(SBS_total_trunk ~ type, random = ~1 | patient, data = df_trunk, method = "ML")
summary(lmm)
```

```
## Linear mixed-effects model fit by maximum likelihood
##   Data: df_trunk 
##        AIC      BIC    logLik
##   801.1926 811.8977 -394.5963
## 
## Random effects:
##  Formula: ~1 | patient
##         (Intercept) Residual
## StdDev:    784.5404 1800.269
## 
## Fixed effects:  SBS_total_trunk ~ type 
##                 Value Std.Error DF   t-value p-value
## (Intercept)  237.3929  804.3788 36 0.2951258  0.7696
## typeadenoma 1247.2940 1099.0622 36 1.1348712  0.2639
## typebc      1534.2804 2082.4552 36 0.7367651  0.4660
## typepolyp   1521.5892  796.1023 36 1.9112986  0.0640
##  Correlation: 
##             (Intr) typdnm typebc
## typeadenoma -0.646              
## typebc      -0.312  0.229       
## typepolyp   -0.781  0.651  0.316
## 
## Standardized Within-Group Residuals:
##        Min         Q1        Med         Q3        Max 
## -1.4853850 -0.4741569 -0.1569824  0.1227567  5.2330544 
## 
## Number of Observations: 44
## Number of Groups: 5
```

``` r
pvalues_lst = c(pvalues_lst, summary(lmm)$tTable[, "p-value"])
```

``` r
lmm <- lme(SBS_total_branch ~ type, random = ~1 | patient/lineage, data = df_trunk_length,
    method = "ML")
summary(lmm)
```

```
## Linear mixed-effects model fit by maximum likelihood
##   Data: df_trunk_length 
##        AIC      BIC    logLik
##   3264.623 3287.816 -1625.312
## 
## Random effects:
##  Formula: ~1 | patient
##         (Intercept)
## StdDev:    148.2895
## 
##  Formula: ~1 | lineage %in% patient
##         (Intercept) Residual
## StdDev:    828.6531 569.2203
## 
## Fixed effects:  SBS_total_branch ~ type 
##                Value Std.Error  DF    t-value p-value
## (Intercept) 813.6191  319.0297 159  2.5502924  0.0117
## typeadenoma -58.9722  462.9690  36 -0.1273782  0.8993
## typebc      413.2742  962.3392  36  0.4294475  0.6702
## typepolyp   975.1467  355.9523  36  2.7395434  0.0095
##  Correlation: 
##             (Intr) typdnm typebc
## typeadenoma -0.684              
## typebc      -0.322  0.222       
## typepolyp   -0.850  0.612  0.288
## 
## Standardized Within-Group Residuals:
##         Min          Q1         Med          Q3         Max 
## -2.75958997 -0.24554079 -0.02500678  0.18418421  6.53762948 
## 
## Number of Observations: 203
## Number of Groups: 
##              patient lineage %in% patient 
##                    5                   44
```

``` r
pvalues_lst = c(pvalues_lst, summary(lmm)$tTable[, "p-value"])
```


``` r
p_corrected = p.adjust(pvalues_lst, method = "BH")
p_table = data.frame(p = pvalues_lst, p_adj = p_corrected)
```
