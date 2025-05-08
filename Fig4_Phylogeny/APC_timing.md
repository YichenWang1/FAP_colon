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
database = read.csv("../data/fap_database_withexposure.csv", header = T)
database$type = as.factor(database$type)
database$type = relevel(database$type, ref = "normal")

# add number of APC mutations
database$num_APC = 3 - (is.na(database$germlinehit) + is.na(database$somatic1) +
    is.na(database$somatic2))
database$num_APC = as.numeric(database$num_APC)

database = database[database$coverage >= 10, ]
database$germlinehit = "APC"

database_exclude_bystander = database[(database$num_APC >= 2 | database$type == "normal"),
    ]
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
##        AIC      BIC    logLik
##   1289.968 1298.069 -641.9838
## 
## Random effects:
##  Formula: ~age - 1 | patient
##              age Residual
## StdDev: 5.070571 66.92812
## 
## Fixed effects:  SBS5 ~ age - 1 
##        Value Std.Error DF  t-value p-value
## age 27.17703  1.357582 14 20.01871       0
## 
## Standardized Within-Group Residuals:
##        Min         Q1        Med         Q3        Max 
## -1.7339257 -0.7818479 -0.0423765  0.4971224  3.9728756 
## 
## Number of Observations: 110
## Number of Groups: 15
```


``` r
df_trunk_length = read.csv("../data/phylogenetic_tree/ACF_polyp_trunk_length.csv")
fixed_slope <- fixef(lmm_sbs5)["age"]
ranef_patient <- ranef(lmm_sbs5) %>%
    tibble::rownames_to_column("patient") %>%
    rename(SBS5_random_slope_patient = age)


# df_trunk_length <- df_trunk_length %>% left_join(ranef_patient, by =
# 'patient') df_trunk_length <- df_trunk_length %>% mutate(SBS5_full_slope =
# fixed_slope + SBS5_random_slope_patient, age_onset_upperlim = SBS5_trunk /
# full_slope) %>% mutate(age_onset_upperlim = if_else(age_onset_upperlim > age,
# age, age_onset_upperlim))

# write.csv(df_trunk_length,'../data/ACF_polyp_trunk_length.csv', quote = F,
# row.names = F)
```


``` r
pvalues_lst = c()
df_trunk_length$trunk_driver_1 = T
df_trunk_length$trunk_driver_1[df_trunk_length$trunk_somatic_driver > 1] = F

lmm <- lme(SBS_total_branch ~ trunk_driver_1, random = ~1 | patient/group, data = df_trunk_length,
    method = "ML")
summary(lmm)  #SBS1
```

```
## Linear mixed-effects model fit by maximum likelihood
##   Data: df_trunk_length 
##        AIC      BIC    logLik
##   1648.612 1661.785 -819.3059
## 
## Random effects:
##  Formula: ~1 | patient
##         (Intercept)
## StdDev:   0.1460727
## 
##  Formula: ~1 | group %in% patient
##         (Intercept) Residual
## StdDev:    968.7454 490.6601
## 
## Fixed effects:  SBS_total_branch ~ trunk_driver_1 
##                        Value Std.Error DF   t-value p-value
## (Intercept)        1787.9826  338.1156 76  5.288080  0.0000
## trunk_driver_1TRUE -381.6785  369.0838 76 -1.034124  0.3044
##  Correlation: 
##                    (Intr)
## trunk_driver_1TRUE -0.807
## 
## Standardized Within-Group Residuals:
##         Min          Q1         Med          Q3         Max 
## -2.23528462 -0.23859288 -0.06358478  0.20587890  6.81741922 
## 
## Number of Observations: 103
## Number of Groups: 
##            patient group %in% patient 
##                  4                 26
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
##        AIC      BIC    logLik
##   513.7046 519.0335 -252.8523
## 
## Random effects:
##  Formula: ~1 | patient
##         (Intercept) Residual
## StdDev:    744.2914 1926.822
## 
## Fixed effects:  SBS_total_trunk ~ trunk_driver_1 
##                        Value Std.Error DF   t-value p-value
## (Intercept)         3788.117  857.3815 23  4.418240  0.0002
## trunk_driver_1TRUE -2984.208  891.8357 23 -3.346141  0.0028
##  Correlation: 
##                    (Intr)
## trunk_driver_1TRUE -0.752
## 
## Standardized Within-Group Residuals:
##        Min         Q1        Med         Q3        Max 
## -1.4883539 -0.2883854 -0.1328292  0.3197319  3.9662225 
## 
## Number of Observations: 28
## Number of Groups: 4
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
##   520.8713 526.2001 -256.4357
## 
## Random effects:
##  Formula: ~1 | patient
##         (Intercept) Residual
## StdDev:    746.8901 2208.127
## 
## Fixed effects:  SBS_total_trunk ~ type 
##                 Value Std.Error DF   t-value p-value
## (Intercept)  533.1182  925.4081 23 0.5760899  0.5701
## typepolyp   1552.5040  970.1615 23 1.6002531  0.1232
##  Correlation: 
##           (Intr)
## typepolyp -0.755
## 
## Standardized Within-Group Residuals:
##        Min         Q1        Med         Q3        Max 
## -1.1888232 -0.5337358 -0.1447241  0.1547360  4.2257192 
## 
## Number of Observations: 28
## Number of Groups: 4
```

``` r
pvalues_lst = c(pvalues_lst, summary(lmm)$tTable[, "p-value"])
```

``` r
lmm <- lme(SBS_total_branch ~ type, random = ~1 | patient/group, data = df_trunk_length,
    method = "ML")
summary(lmm)
```

```
## Linear mixed-effects model fit by maximum likelihood
##   Data: df_trunk_length 
##        AIC    BIC    logLik
##   1645.127 1658.3 -817.5633
## 
## Random effects:
##  Formula: ~1 | patient
##         (Intercept)
## StdDev:   0.2648246
## 
##  Formula: ~1 | group %in% patient
##         (Intercept) Residual
## StdDev:    946.4546 483.3642
## 
## Fixed effects:  SBS_total_branch ~ type 
##                Value Std.Error DF  t-value p-value
## (Intercept) 878.8537  350.5046 77 2.507396  0.0143
## typepolyp   908.1280  421.8047 21 2.152959  0.0431
##  Correlation: 
##           (Intr)
## typepolyp -0.831
## 
## Standardized Within-Group Residuals:
##         Min          Q1         Med          Q3         Max 
## -2.26556920 -0.22622869 -0.03054915  0.17370788  6.92281500 
## 
## Number of Observations: 103
## Number of Groups: 
##            patient group %in% patient 
##                  4                 26
```

``` r
pvalues_lst = c(pvalues_lst, summary(lmm)$tTable[, "p-value"])
```


``` r
p_corrected = p.adjust(pvalues_lst, method = "BH")
p_table = data.frame(p = pvalues_lst, p_adj = p_corrected)
```
