---
title: "hypothesis_testing"
author: "Yichen Wang"
output:   
  html_document:
    keep_md: yes
---




``` r
library(ggplot2)
library(RColorBrewer)
library(nlme)
library(tidyverse)
```


``` r
dot_cols = c("#DF5B3F", "#F5B46F", "#70A3C4", "black")
names(dot_cols) = c("Small and large polyp", "Diminutive polyp", "ACF", "Age effect")


df_sbs_disease = data.frame(matrix(nrow = 12, ncol = 5))
colnames(df_sbs_disease) = c("rate", "type", "group", "lower_95_CI", "upper_95_CI")
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




# Whether germline APC mutation leads to incresed sbs burden


``` r
df = full_dataset[full_dataset$type %in% c("normal"), ]
lmm_sbs <- lme(sbstotal_corr ~ age:germlinehit, random = ~1 | patient/block, data = df,
    method = "ML")
lmm2_sbs <- lme(sbstotal_corr ~ age, random = ~1 | patient/block, data = df, method = "ML")
anova(lmm_sbs, lmm2_sbs)
```

```
##          Model df      AIC      BIC    logLik   Test   L.Ratio p-value
## lmm_sbs      1  6 7669.979 7695.397 -3828.990                         
## lmm2_sbs     2  5 7668.273 7689.455 -3829.137 1 vs 2 0.2943922  0.5874
```


``` r
df = database_exclude_bystander[database_exclude_bystander$type %in% c("normal"),
    ]
lmm2_sbs <- lme(sbstotal_corr ~ age, random = ~1 | patient/block, data = df, method = "ML")
summary(lmm2_sbs)
```

```
## Linear mixed-effects model fit by maximum likelihood
##   Data: df 
##        AIC      BIC    logLik
##   1487.407 1500.909 -738.7034
## 
## Random effects:
##  Formula: ~1 | patient
##         (Intercept)
## StdDev:    273.9314
## 
##  Formula: ~1 | block %in% patient
##         (Intercept) Residual
## StdDev:   0.0217301 166.6321
## 
## Fixed effects:  sbstotal_corr ~ age 
##                Value Std.Error DF  t-value p-value
## (Intercept) 391.4692 267.28280 88 1.464626  0.1466
## age          41.7334   8.77144 13 4.757877  0.0004
##  Correlation: 
##     (Intr)
## age -0.96 
## 
## Standardized Within-Group Residuals:
##        Min         Q1        Med         Q3        Max 
## -1.9792678 -0.5219872 -0.1569062  0.4917597  3.0142858 
## 
## Number of Observations: 110
## Number of Groups: 
##            patient block %in% patient 
##                 15                 22
```

``` r
intervals(lmm2_sbs, which = "fixed")
```

```
## Approximate 95% confidence intervals
## 
##  Fixed effects:
##                  lower      est.    upper
## (Intercept) -134.84822 391.46921 917.7866
## age           22.95694  41.73342  60.5099
```

``` r
pvalues_lst = c(pvalues_lst, summary(lmm2_sbs)$tTable[, "p-value"])
```


# Whether germline APC mutation leads to incresed Indel burden


``` r
df = full_dataset[full_dataset$type %in% c("normal"), ]
lmm_id <- lme(indeltotal_corr ~ age:germlinehit, random = ~1 | patient/block, data = df,
    method = "ML")
lmm2_id <- lme(indeltotal_corr ~ age, random = ~1 | patient/block, data = df, method = "ML")
anova(lmm_id, lmm2_id)
```

```
##         Model df     AIC      BIC   logLik   Test  L.Ratio p-value
## lmm_id      1  6 4988.12 5013.538 -2488.06                        
## lmm2_id     2  5 4987.64 5008.822 -2488.82 1 vs 2 1.520033  0.2176
```



``` r
df = database_exclude_bystander[database_exclude_bystander$type %in% c("normal"),
    ]
lmm2_id <- lme(indeltotal_corr ~ age, random = ~1 | patient, data = df, method = "ML")
summary(lmm2_id)
```

```
## Linear mixed-effects model fit by maximum likelihood
##   Data: df 
##        AIC      BIC   logLik
##   949.6499 960.4518 -470.825
## 
## Random effects:
##  Formula: ~1 | patient
##         (Intercept) Residual
## StdDev:    29.81262 14.13942
## 
## Fixed effects:  indeltotal_corr ~ age 
##                Value Std.Error DF  t-value p-value
## (Intercept) 45.20409 28.658991 95 1.577309  0.1180
## age          1.74830  0.940567 13 1.858774  0.0858
##  Correlation: 
##     (Intr)
## age -0.96 
## 
## Standardized Within-Group Residuals:
##         Min          Q1         Med          Q3         Max 
## -2.63825904 -0.57898587 -0.04102255  0.56023221  2.10693577 
## 
## Number of Observations: 110
## Number of Groups: 15
```

``` r
fixed.m1 <- data.frame(fixef(lmm2_id))
intervals(lmm2_id, which = "fixed")
```

```
## Approximate 95% confidence intervals
## 
##  Fixed effects:
##                   lower      est.      upper
## (Intercept) -11.1715950 45.204093 101.579782
## age          -0.2651129  1.748301   3.761714
```

``` r
pvalues_lst = c(pvalues_lst, summary(lmm2_id)$tTable[, "p-value"])
```


``` r
p_corrected = p.adjust(pvalues_lst, method = "BH")
p_table = data.frame(p = pvalues_lst, p_adj = p_corrected)
```

## Whetehr bystander crypts are the same as normal




``` r
df = full_dataset_bystander
lmm_sbs <- lme(sbstotal_corr ~ age + germlinehit, random = ~1 | patient/block/group_id,
    data = df, method = "ML")
summary(lmm_sbs)
```

```
## Linear mixed-effects model fit by maximum likelihood
##   Data: df 
##        AIC      BIC    logLik
##   8472.012 8502.432 -4229.006
## 
## Random effects:
##  Formula: ~1 | patient
##         (Intercept)
## StdDev:    433.3903
## 
##  Formula: ~1 | block %in% patient
##         (Intercept)
## StdDev:    367.1388
## 
##  Formula: ~1 | group_id %in% block %in% patient
##         (Intercept) Residual
## StdDev:    303.6036 143.7448
## 
## Fixed effects:  sbstotal_corr ~ age + germlinehit 
##                   Value Std.Error  DF  t-value p-value
## (Intercept)    418.0751  350.2753 410 1.193562  0.2333
## age             35.9782    5.5940  53 6.431548  0.0000
## germlinehitAPC 104.1120  242.8081  53 0.428783  0.6698
##  Correlation: 
##                (Intr) age   
## age            -0.972       
## germlinehitAPC -0.787  0.728
## 
## Standardized Within-Group Residuals:
##         Min          Q1         Med          Q3         Max 
## -1.82525656 -0.20666223 -0.01525621  0.20012953  2.23493878 
## 
## Number of Observations: 570
## Number of Groups: 
##                          patient               block %in% patient 
##                               56                              123 
## group_id %in% block %in% patient 
##                              533
```


``` r
df = full_dataset_bystander
lmm_id <- lme(indeltotal_corr ~ age + germlinehit, random = ~1 | patient/block/group_id,
    data = df, method = "ML")
summary(lmm_id)
```

```
## Linear mixed-effects model fit by maximum likelihood
##   Data: df 
##        AIC      BIC    logLik
##   5502.201 5532.621 -2744.101
## 
## Random effects:
##  Formula: ~1 | patient
##         (Intercept)
## StdDev:     31.8553
## 
##  Formula: ~1 | block %in% patient
##         (Intercept)
## StdDev:    21.98671
## 
##  Formula: ~1 | group_id %in% block %in% patient
##         (Intercept) Residual
## StdDev:    21.04132  13.6108
## 
## Fixed effects:  indeltotal_corr ~ age + germlinehit 
##                   Value Std.Error  DF  t-value p-value
## (Intercept)    47.93071 24.295841 410 1.972795  0.0492
## age             0.96404  0.388766  53 2.479732  0.0164
## germlinehitAPC 18.74709 16.889145  53 1.110008  0.2720
##  Correlation: 
##                (Intr) age   
## age            -0.971       
## germlinehitAPC -0.784  0.723
## 
## Standardized Within-Group Residuals:
##        Min         Q1        Med         Q3        Max 
## -2.4409627 -0.2957849 -0.0341705  0.2553690  2.9022896 
## 
## Number of Observations: 570
## Number of Groups: 
##                          patient               block %in% patient 
##                               56                              123 
## group_id %in% block %in% patient 
##                              533
```


# Whether ACF leads to incresed sbs burden

``` r
pvalues_lst = c()
df = full_dataset[full_dataset$type %in% c("normal", "acf"), ]
lmm_sbs18 <- lme(SBS18 ~ age + type, random = ~1 | patient/block/group_id, data = df,
    method = "ML")
summary(lmm_sbs18)
```

```
## Linear mixed-effects model fit by maximum likelihood
##   Data: df 
##        AIC      BIC    logLik
##   6768.002 6798.172 -3377.001
## 
## Random effects:
##  Formula: ~1 | patient
##         (Intercept)
## StdDev:    109.9643
## 
##  Formula: ~1 | block %in% patient
##         (Intercept)
## StdDev:    113.3981
## 
##  Formula: ~1 | group_id %in% block %in% patient
##         (Intercept) Residual
## StdDev:    88.88031 29.56008
## 
## Fixed effects:  SBS18 ~ age + type 
##                 Value Std.Error  DF  t-value p-value
## (Intercept)   3.23063  59.80029 402 0.054024  0.9569
## age           5.82504   1.05354  54 5.528994  0.0000
## typeacf     127.86976  45.82769 402 2.790229  0.0055
##  Correlation: 
##         (Intr) age   
## age     -0.946       
## typeacf -0.055  0.042
## 
## Standardized Within-Group Residuals:
##          Min           Q1          Med           Q3          Max 
## -2.392746558 -0.152123371 -0.007844915  0.129410192  2.746106871 
## 
## Number of Observations: 550
## Number of Groups: 
##                          patient               block %in% patient 
##                               56                              118 
## group_id %in% block %in% patient 
##                              521
```

``` r
pvalues_lst = c(pvalues_lst, summary(lmm_sbs18)$tTable[, "p-value"])
intervals(lmm_sbs18, which = "fixed")
```

```
## Approximate 95% confidence intervals
## 
##  Fixed effects:
##                   lower       est.      upper
## (Intercept) -114.008660   3.230630 120.469921
## age            3.718581   5.825042   7.931503
## typeacf       38.023941 127.869755 217.715570
```



``` r
df = full_dataset[full_dataset$type %in% c("normal", "acf"), ]
lmm_sbs1 <- lme(SBS1 ~ age + type, random = ~1 | patient/block/group_id, data = df,
    method = "ML")
summary(lmm_sbs1)
```

```
## Linear mixed-effects model fit by maximum likelihood
##   Data: df 
##       AIC      BIC    logLik
##   7359.49 7389.659 -3672.745
## 
## Random effects:
##  Formula: ~1 | patient
##         (Intercept)
## StdDev:    165.6206
## 
##  Formula: ~1 | block %in% patient
##         (Intercept)
## StdDev:    172.2647
## 
##  Formula: ~1 | group_id %in% block %in% patient
##         (Intercept) Residual
## StdDev:    153.9326 56.09518
## 
## Fixed effects:  SBS1 ~ age + type 
##                 Value Std.Error  DF  t-value p-value
## (Intercept) 163.18929  91.68268 402 1.779936  0.0758
## age          11.84235   1.61463  54 7.334384  0.0000
## typeacf     152.75086  77.63271 402 1.967609  0.0498
##  Correlation: 
##         (Intr) age   
## age     -0.946       
## typeacf -0.060  0.046
## 
## Standardized Within-Group Residuals:
##          Min           Q1          Med           Q3          Max 
## -1.896643930 -0.186525097 -0.007174227  0.162050951  2.492872614 
## 
## Number of Observations: 550
## Number of Groups: 
##                          patient               block %in% patient 
##                               56                              118 
## group_id %in% block %in% patient 
##                              521
```

``` r
pvalues_lst = c(pvalues_lst, summary(lmm_sbs1)$tTable[, "p-value"])
intervals(lmm_sbs1, which = "fixed")
```

```
## Approximate 95% confidence intervals
## 
##  Fixed effects:
##                   lower      est.     upper
## (Intercept) -16.5558753 163.18929 342.93445
## age           8.6140432  11.84235  15.07066
## typeacf       0.5508541 152.75086 304.95086
```


``` r
df = full_dataset[full_dataset$type %in% c("normal", "acf"), ]
lmm_sbs5 <- lme(SBS5 ~ age + type, random = ~1 | patient/block/group_id, data = df,
    method = "ML")
summary(lmm_sbs5)
```

```
## Linear mixed-effects model fit by maximum likelihood
##   Data: df 
##        AIC      BIC    logLik
##   7066.953 7097.123 -3526.477
## 
## Random effects:
##  Formula: ~1 | patient
##         (Intercept)
## StdDev:    167.9764
## 
##  Formula: ~1 | block %in% patient
##         (Intercept)
## StdDev:    135.5575
## 
##  Formula: ~1 | group_id %in% block %in% patient
##         (Intercept) Residual
## StdDev:    113.4226 44.72813
## 
## Fixed effects:  SBS5 ~ age + type 
##                 Value Std.Error  DF   t-value p-value
## (Intercept) 238.18854  82.81507 402  2.876150  0.0042
## age          16.29138   1.46785  54 11.098803  0.0000
## typeacf      25.03663  58.59092 402  0.427312  0.6694
##  Correlation: 
##         (Intr) age   
## age     -0.944       
## typeacf -0.046  0.035
## 
## Standardized Within-Group Residuals:
##          Min           Q1          Med           Q3          Max 
## -1.691210446 -0.186299016 -0.005027297  0.179011226  2.706483229 
## 
## Number of Observations: 550
## Number of Groups: 
##                          patient               block %in% patient 
##                               56                              118 
## group_id %in% block %in% patient 
##                              521
```

``` r
pvalues_lst = c(pvalues_lst, summary(lmm_sbs5)$tTable[, "p-value"])
intervals(lmm_sbs5, which = "fixed")
```

```
## Approximate 95% confidence intervals
## 
##  Fixed effects:
##                 lower      est.    upper
## (Intercept)  75.82846 238.18854 400.5486
## age          13.35655  16.29138  19.2262
## typeacf     -89.83169  25.03663 139.9049
```


``` r
df = full_dataset[full_dataset$type %in% c("normal", "acf"), ]
lmm_sbs <- lme(sbstotal_corr ~ age + type, random = ~1 | patient/block/group_id,
    data = df, method = "ML")
summary(lmm_sbs)
```

```
## Linear mixed-effects model fit by maximum likelihood
##   Data: df 
##        AIC      BIC    logLik
##   8184.295 8214.465 -4085.148
## 
## Random effects:
##  Formula: ~1 | patient
##         (Intercept)
## StdDev:    421.8872
## 
##  Formula: ~1 | block %in% patient
##         (Intercept)
## StdDev:    386.0577
## 
##  Formula: ~1 | group_id %in% block %in% patient
##         (Intercept) Residual
## StdDev:    322.7354 109.2594
## 
## Fixed effects:  sbstotal_corr ~ age + type 
##                Value Std.Error  DF  t-value p-value
## (Intercept) 567.9993 218.52879 402 2.599196  0.0097
## age          33.7659   3.86121  54 8.744905  0.0000
## typeacf     337.9887 164.96785 402 2.048816  0.0411
##  Correlation: 
##         (Intr) age   
## age     -0.945       
## typeacf -0.051  0.039
## 
## Standardized Within-Group Residuals:
##         Min          Q1         Med          Q3         Max 
## -1.34711080 -0.16700088 -0.01247547  0.14861225  2.97073714 
## 
## Number of Observations: 550
## Number of Groups: 
##                          patient               block %in% patient 
##                               56                              118 
## group_id %in% block %in% patient 
##                              521
```

``` r
pvalues_lst = c(pvalues_lst, summary(lmm_sbs)$tTable[, "p-value"])
intervals(lmm_sbs, which = "fixed")
```

```
## Approximate 95% confidence intervals
## 
##  Fixed effects:
##                 lower      est.     upper
## (Intercept) 139.57055 567.99927 996.42799
## age          26.04578  33.76589  41.48601
## typeacf      14.56695 337.98871 661.41046
```

# Whether ACF leads to incresed Indel burden

``` r
df = full_dataset[full_dataset$type %in% c("normal", "acf"), ]
lmm_id1 <- lme(ID1 ~ age + type, random = ~1 | patient/block/group_id, data = df,
    method = "ML")
summary(lmm_id1)
```

```
## Linear mixed-effects model fit by maximum likelihood
##   Data: df 
##        AIC      BIC    logLik
##   4244.349 4274.518 -2115.174
## 
## Random effects:
##  Formula: ~1 | patient
##         (Intercept)
## StdDev:    8.814155
## 
##  Formula: ~1 | block %in% patient
##         (Intercept)
## StdDev:    7.654835
## 
##  Formula: ~1 | group_id %in% block %in% patient
##         (Intercept) Residual
## StdDev:    7.370831 6.298343
## 
## Fixed effects:  ID1 ~ age + type 
##                 Value Std.Error  DF  t-value p-value
## (Intercept) 11.964024  4.649177 402 2.573364  0.0104
## age          0.286909  0.082179  54 3.491257  0.0010
## typeacf     28.789365  4.049845 402 7.108758  0.0000
##  Correlation: 
##         (Intr) age   
## age     -0.945       
## typeacf -0.058  0.044
## 
## Standardized Within-Group Residuals:
##         Min          Q1         Med          Q3         Max 
## -4.31626477 -0.31530900 -0.05395308  0.29082506  2.94121241 
## 
## Number of Observations: 550
## Number of Groups: 
##                          patient               block %in% patient 
##                               56                              118 
## group_id %in% block %in% patient 
##                              521
```

``` r
pvalues_lst = c(pvalues_lst, summary(lmm_id1)$tTable[, "p-value"])
intervals(lmm_id1, which = "fixed")
```

```
## Approximate 95% confidence intervals
## 
##  Fixed effects:
##                  lower      est.      upper
## (Intercept)  2.8492482 11.964024 21.0788001
## age          0.1225994  0.286909  0.4512186
## typeacf     20.8495885 28.789365 36.7291412
```


``` r
df = full_dataset[full_dataset$type %in% c("normal", "acf"), ]
lmm_id2 <- lme(ID2 ~ age + type, random = ~1 | patient/block/group_id, data = df,
    method = "ML")
summary(lmm_id2)
```

```
## Linear mixed-effects model fit by maximum likelihood
##   Data: df 
##        AIC      BIC    logLik
##   3674.931 3705.101 -1830.466
## 
## Random effects:
##  Formula: ~1 | patient
##         (Intercept)
## StdDev:    4.190909
## 
##  Formula: ~1 | block %in% patient
##         (Intercept)
## StdDev:     5.05909
## 
##  Formula: ~1 | group_id %in% block %in% patient
##         (Intercept) Residual
## StdDev:    5.676997 2.012016
## 
## Fixed effects:  ID2 ~ age + type 
##                 Value Std.Error  DF  t-value p-value
## (Intercept)  6.536366 2.5915675 402 2.522167  0.0120
## age          0.115753 0.0454909  54 2.544537  0.0138
## typeacf     11.040299 2.7041790 402 4.082680  0.0001
##  Correlation: 
##         (Intr) age   
## age     -0.948       
## typeacf -0.074  0.057
## 
## Standardized Within-Group Residuals:
##         Min          Q1         Med          Q3         Max 
## -1.86943951 -0.17077905 -0.01231143  0.16702702  1.74428706 
## 
## Number of Observations: 550
## Number of Groups: 
##                          patient               block %in% patient 
##                               56                              118 
## group_id %in% block %in% patient 
##                              521
```

``` r
pvalues_lst = c(pvalues_lst, summary(lmm_id2)$tTable[, "p-value"])
intervals(lmm_id2, which = "fixed")
```

```
## Approximate 95% confidence intervals
## 
##  Fixed effects:
##                 lower       est.      upper
## (Intercept) 1.4555621  6.5363661 11.6171701
## age         0.0247986  0.1157532  0.2067079
## typeacf     5.7387185 11.0402988 16.3418791
```


``` r
df = full_dataset[full_dataset$type %in% c("normal", "acf"), ]
lmm_id5 <- lme(ID5 ~ age + type, random = ~1 | patient/block/group_id, data = df,
    method = "ML")
summary(lmm_id5)
```

```
## Linear mixed-effects model fit by maximum likelihood
##   Data: df 
##        AIC      BIC    logLik
##   4432.373 4462.542 -2209.186
## 
## Random effects:
##  Formula: ~1 | patient
##         (Intercept)
## StdDev:    13.26823
## 
##  Formula: ~1 | block %in% patient
##         (Intercept)
## StdDev:    5.686359
## 
##  Formula: ~1 | group_id %in% block %in% patient
##         (Intercept) Residual
## StdDev:    11.32731 4.368708
## 
## Fixed effects:  ID5 ~ age + type 
##                 Value Std.Error  DF  t-value p-value
## (Intercept) 22.804261  5.942143 402 3.837717  0.0001
## age          0.437612  0.106241  54 4.119071  0.0001
## typeacf      3.321770  4.888164 402 0.679554  0.4972
##  Correlation: 
##         (Intr) age   
## age     -0.942       
## typeacf -0.037  0.027
## 
## Standardized Within-Group Residuals:
##           Min            Q1           Med            Q3           Max 
## -2.5042391618 -0.2048001779  0.0008482901  0.2007231895  2.6117672427 
## 
## Number of Observations: 550
## Number of Groups: 
##                          patient               block %in% patient 
##                               56                              118 
## group_id %in% block %in% patient 
##                              521
```

``` r
pvalues_lst = c(pvalues_lst, summary(lmm_id5)$tTable[, "p-value"])
intervals(lmm_id5, which = "fixed")
```

```
## Approximate 95% confidence intervals
## 
##  Fixed effects:
##                  lower       est.      upper
## (Intercept) 11.1546071 22.8042611 34.4539150
## age          0.2251946  0.4376124  0.6500302
## typeacf     -6.2615438  3.3217697 12.9050831
```


``` r
df = full_dataset[full_dataset$type %in% c("normal", "acf"), ]
lmm_id <- lme(indeltotal_corr ~ age + type, random = ~1 | patient/block/group_id,
    data = df, method = "ML")
summary(lmm_id)
```

```
## Linear mixed-effects model fit by maximum likelihood
##   Data: df 
##        AIC      BIC    logLik
##   5304.629 5334.799 -2645.315
## 
## Random effects:
##  Formula: ~1 | patient
##         (Intercept)
## StdDev:    31.26036
## 
##  Formula: ~1 | block %in% patient
##         (Intercept)
## StdDev:    23.98872
## 
##  Formula: ~1 | group_id %in% block %in% patient
##         (Intercept) Residual
## StdDev:    23.90356 8.336887
## 
## Fixed effects:  indeltotal_corr ~ age + type 
##                Value Std.Error  DF  t-value p-value
## (Intercept) 70.61719 15.335412 402 4.604845  0.0000
## age          0.62959  0.272018  54 2.314517  0.0245
## typeacf     35.13579 11.850217 402 2.964991  0.0032
##  Correlation: 
##         (Intr) age   
## age     -0.944       
## typeacf -0.048  0.036
## 
## Standardized Within-Group Residuals:
##         Min          Q1         Med          Q3         Max 
## -3.07908973 -0.16943069 -0.01480108  0.15752801  2.20585344 
## 
## Number of Observations: 550
## Number of Groups: 
##                          patient               block %in% patient 
##                               56                              118 
## group_id %in% block %in% patient 
##                              521
```

``` r
pvalues_lst = c(pvalues_lst, summary(lmm_id)$tTable[, "p-value"])
intervals(lmm_id, which = "fixed")
```

```
## Approximate 95% confidence intervals
## 
##  Fixed effects:
##                   lower       est.      upper
## (Intercept) 40.55190774 70.6171943 100.682481
## age          0.08571639  0.6295909   1.173465
## typeacf     11.90327796 35.1357920  58.368306
```

``` r
p_corrected_acf = p.adjust(pvalues_lst, method = "BH")
p_table_acf = data.frame(p = pvalues_lst, p_adj = p_corrected_acf)
```

# Whether Polyp leads to incresed sbs burden

``` r
pvalues_lst = c()
df = full_dataset[full_dataset$type %in% c("normal", "polyp", "acf"), ]
lmm_sbs18 <- lme(SBS18 ~ type + age, random = ~1 | patient/block/group_id, data = df,
    method = "ML")
summary(lmm_sbs18)
```

```
## Linear mixed-effects model fit by maximum likelihood
##   Data: df 
##        AIC      BIC    logLik
##   7586.593 7621.728 -3785.297
## 
## Random effects:
##  Formula: ~1 | patient
##         (Intercept)
## StdDev:    111.1538
## 
##  Formula: ~1 | block %in% patient
##         (Intercept)
## StdDev:    109.0696
## 
##  Formula: ~1 | group_id %in% block %in% patient
##         (Intercept) Residual
## StdDev:    105.1554 60.49361
## 
## Fixed effects:  SBS18 ~ type + age 
##                Value Std.Error  DF   t-value p-value
## (Intercept)  -6.8096  60.63126 409 -0.112312  0.9106
## typeacf      82.0658  49.10143 409  1.671352  0.0954
## typepolyp   551.8790  46.94625 409 11.755550  0.0000
## age           5.9829   1.07031  54  5.589837  0.0000
##  Correlation: 
##           (Intr) typecf typply
## typeacf   -0.082              
## typepolyp -0.134  0.356       
## age       -0.946  0.066  0.113
## 
## Standardized Within-Group Residuals:
##         Min          Q1         Med          Q3         Max 
## -2.94886898 -0.20233661 -0.01166645  0.17841854  4.86175978 
## 
## Number of Observations: 597
## Number of Groups: 
##                          patient               block %in% patient 
##                               56                              125 
## group_id %in% block %in% patient 
##                              536
```

``` r
pvalues_lst = c(pvalues_lst, summary(lmm_sbs18)$tTable[, "p-value"])
intervals(lmm_sbs18, which = "fixed")
```

```
## Approximate 95% confidence intervals
## 
##  Fixed effects:
##                   lower       est.      upper
## (Intercept) -125.597408  -6.809593 111.978222
## typeacf      -14.132969  82.065795 178.264559
## typepolyp    459.902607 551.878970 643.855333
## age            3.844219   5.982862   8.121505
```



``` r
df = full_dataset[full_dataset$type %in% c("normal", "polyp", "acf"), ]
lmm_sbs1 <- lme(SBS1 ~ age + type, random = ~1 | patient/block/group_id, data = df,
    method = "ML")
summary(lmm_sbs1)
```

```
## Linear mixed-effects model fit by maximum likelihood
##   Data: df 
##        AIC      BIC    logLik
##   8002.436 8037.572 -3993.218
## 
## Random effects:
##  Formula: ~1 | patient
##         (Intercept)
## StdDev:    170.5543
## 
##  Formula: ~1 | block %in% patient
##         (Intercept)
## StdDev:    166.0696
## 
##  Formula: ~1 | group_id %in% block %in% patient
##         (Intercept) Residual
## StdDev:     158.044 71.11303
## 
## Fixed effects:  SBS1 ~ age + type 
##                Value Std.Error  DF  t-value p-value
## (Intercept) 150.1547  91.94321 409 1.633125  0.1032
## age          12.0610   1.62356  54 7.428699  0.0000
## typeacf     108.2394  72.30196 409 1.497047  0.1352
## typepolyp   562.7483  69.65983 409 8.078519  0.0000
##  Correlation: 
##           (Intr) age    typecf
## age       -0.945              
## typeacf   -0.079  0.064       
## typepolyp -0.129  0.109  0.358
## 
## Standardized Within-Group Residuals:
##          Min           Q1          Med           Q3          Max 
## -2.313293592 -0.214895088 -0.007314206  0.199859935  2.876333599 
## 
## Number of Observations: 597
## Number of Groups: 
##                          patient               block %in% patient 
##                               56                              125 
## group_id %in% block %in% patient 
##                              536
```

``` r
pvalues_lst = c(pvalues_lst, summary(lmm_sbs1)$tTable[, "p-value"])
intervals(lmm_sbs1, which = "fixed")
```

```
## Approximate 95% confidence intervals
## 
##  Fixed effects:
##                  lower      est.    upper
## (Intercept) -29.978969 150.15473 330.2884
## age           8.816847  12.06098  15.3051
## typeacf     -33.413462 108.23941 249.8923
## typepolyp   426.271804 562.74827 699.2247
```



``` r
df = full_dataset[full_dataset$type %in% c("normal", "polyp", "acf"), ]
lmm_sbs5 <- lme(SBS5 ~ age + type, random = ~1 | patient/block/group_id, data = df,
    method = "ML")
summary(lmm_sbs5)
```

```
## Linear mixed-effects model fit by maximum likelihood
##   Data: df 
##        AIC      BIC    logLik
##   7887.831 7922.966 -3935.915
## 
## Random effects:
##  Formula: ~1 | patient
##         (Intercept)
## StdDev:    165.0133
## 
##  Formula: ~1 | block %in% patient
##         (Intercept)
## StdDev:    116.4379
## 
##  Formula: ~1 | group_id %in% block %in% patient
##         (Intercept) Residual
## StdDev:    152.4777 59.51892
## 
## Fixed effects:  SBS5 ~ age + type 
##                Value Std.Error  DF   t-value p-value
## (Intercept) 227.7255  81.00279 409  2.811329  0.0052
## age          16.4383   1.43870  54 11.425778  0.0000
## typeacf     -81.8915  66.18808 409 -1.237254  0.2167
## typepolyp   398.6964  60.76248 409  6.561555  0.0000
##  Correlation: 
##           (Intr) age    typecf
## age       -0.944              
## typeacf   -0.063  0.050       
## typepolyp -0.103  0.086  0.326
## 
## Standardized Within-Group Residuals:
##          Min           Q1          Med           Q3          Max 
## -2.481238294 -0.142773971 -0.009913201  0.155023376  3.698629077 
## 
## Number of Observations: 597
## Number of Groups: 
##                          patient               block %in% patient 
##                               56                              125 
## group_id %in% block %in% patient 
##                              536
```

``` r
pvalues_lst = c(pvalues_lst, summary(lmm_sbs5)$tTable[, "p-value"])
intervals(lmm_sbs5, which = "fixed")
```

```
## Approximate 95% confidence intervals
## 
##  Fixed effects:
##                  lower      est.     upper
## (Intercept)   69.02611 227.72552 386.42494
## age           13.56356  16.43831  19.31306
## typeacf     -211.56612 -81.89146  47.78321
## typepolyp    279.65145 398.69636 517.74127
```



``` r
df = full_dataset[full_dataset$type %in% c("normal", "polyp", "acf"), ]
lmm_sbs <- lme(sbstotal_corr ~ age + type, random = ~1 | patient/block/group_id,
    data = df, method = "ML")
summary(lmm_sbs)
```

```
## Linear mixed-effects model fit by maximum likelihood
##   Data: df 
##        AIC      BIC    logLik
##   9135.383 9170.519 -4559.692
## 
## Random effects:
##  Formula: ~1 | patient
##         (Intercept)
## StdDev:    422.0505
## 
##  Formula: ~1 | block %in% patient
##         (Intercept)
## StdDev:    344.5139
## 
##  Formula: ~1 | group_id %in% block %in% patient
##         (Intercept) Residual
## StdDev:    439.9394 163.5798
## 
## Fixed effects:  sbstotal_corr ~ age + type 
##                 Value Std.Error  DF  t-value p-value
## (Intercept)  531.7227 218.48084 409 2.433727  0.0154
## age           34.3347   3.86886  54 8.874646  0.0000
## typeacf       26.3095 190.59467 409 0.138039  0.8903
## typepolyp   1580.0229 175.46502 409 9.004774  0.0000
##  Correlation: 
##           (Intr) age    typecf
## age       -0.945              
## typeacf   -0.072  0.057       
## typepolyp -0.118  0.099  0.324
## 
## Standardized Within-Group Residuals:
##         Min          Q1         Med          Q3         Max 
## -2.80194466 -0.14885323 -0.00535994  0.15037503  3.87921095 
## 
## Number of Observations: 597
## Number of Groups: 
##                          patient               block %in% patient 
##                               56                              125 
## group_id %in% block %in% patient 
##                              536
```

``` r
pvalues_lst = c(pvalues_lst, summary(lmm_sbs)$tTable[, "p-value"])
intervals(lmm_sbs, which = "fixed")
```

```
## Approximate 95% confidence intervals
## 
##  Fixed effects:
##                  lower       est.     upper
## (Intercept)  103.67846  531.72274  959.7670
## age           26.60417   34.33474   42.0653
## typeacf     -347.10065   26.30947  399.7196
## typepolyp   1236.25456 1580.02290 1923.7912
```


# Whether polyps leads to incresed Indel burden

``` r
df = full_dataset[full_dataset$type %in% c("normal", "polyp", "acf"), ]
lmm_id1 <- lme(ID1 ~ age + type, random = ~1 | patient/block/group_id, data = df,
    method = "ML")
summary(lmm_id1)
```

```
## Linear mixed-effects model fit by maximum likelihood
##   Data: df 
##        AIC      BIC    logLik
##   5155.947 5191.083 -2569.974
## 
## Random effects:
##  Formula: ~1 | patient
##         (Intercept)
## StdDev:    11.54791
## 
##  Formula: ~1 | block %in% patient
##         (Intercept)
## StdDev:    12.43906
## 
##  Formula: ~1 | group_id %in% block %in% patient
##         (Intercept) Residual
## StdDev:    13.75879 8.580074
## 
## Fixed effects:  ID1 ~ age + type 
##                Value Std.Error  DF   t-value p-value
## (Intercept) 11.30530  6.772683 409  1.669250  0.0958
## age          0.30151  0.119286  54  2.527650  0.0144
## typeacf     24.09452  6.327126 409  3.808131  0.0002
## typepolyp   96.54454  5.884710 409 16.405998  0.0000
##  Correlation: 
##           (Intr) age    typecf
## age       -0.947              
## typeacf   -0.093  0.075       
## typepolyp -0.153  0.129  0.334
## 
## Standardized Within-Group Residuals:
##        Min         Q1        Med         Q3        Max 
## -3.7095248 -0.1777413 -0.0153875  0.1723019  4.1690180 
## 
## Number of Observations: 597
## Number of Groups: 
##                          patient               block %in% patient 
##                               56                              125 
## group_id %in% block %in% patient 
##                              536
```

``` r
pvalues_lst = c(pvalues_lst, summary(lmm_id1)$tTable[, "p-value"])
intervals(lmm_id1, which = "fixed")
```

```
## Approximate 95% confidence intervals
## 
##  Fixed effects:
##                   lower      est.       upper
## (Intercept) -1.96363695 11.305299  24.5742355
## age          0.06316167  0.301513   0.5398643
## typeacf     11.69851606 24.094522  36.4905285
## typepolyp   85.01531054 96.544543 108.0737755
```



``` r
df = full_dataset[full_dataset$type %in% c("normal", "polyp", "acf"), ]
lmm_id2 <- lme(ID2 ~ age + type, random = ~1 | patient/block/group_id, data = df,
    method = "ML")
summary(lmm_id2)
```

```
## Linear mixed-effects model fit by maximum likelihood
##   Data: df 
##        AIC      BIC    logLik
##   4881.567 4916.702 -2432.783
## 
## Random effects:
##  Formula: ~1 | patient
##         (Intercept)
## StdDev:    4.088115
## 
##  Formula: ~1 | block %in% patient
##         (Intercept)
## StdDev:    2.379808
## 
##  Formula: ~1 | group_id %in% block %in% patient
##         (Intercept) Residual
## StdDev:     13.3938 6.262595
## 
## Fixed effects:  ID2 ~ age + type 
##                Value Std.Error  DF   t-value p-value
## (Intercept)  4.46680  3.027102 409  1.475603  0.1408
## age          0.15244  0.053862  54  2.830126  0.0065
## typeacf      7.04493  4.906475 409  1.435844  0.1518
## typepolyp   51.75983  4.040676 409 12.809698  0.0000
##  Correlation: 
##           (Intr) age    typecf
## age       -0.947              
## typeacf   -0.101  0.078       
## typepolyp -0.163  0.135  0.146
## 
## Standardized Within-Group Residuals:
##         Min          Q1         Med          Q3         Max 
## -5.55355768 -0.13855274 -0.02077835  0.11597674  4.06075793 
## 
## Number of Observations: 597
## Number of Groups: 
##                          patient               block %in% patient 
##                               56                              125 
## group_id %in% block %in% patient 
##                              536
```

``` r
pvalues_lst = c(pvalues_lst, summary(lmm_id2)$tTable[, "p-value"])
intervals(lmm_id2, which = "fixed")
```

```
## Approximate 95% confidence intervals
## 
##  Fixed effects:
##                   lower       est.      upper
## (Intercept) -1.46384917  4.4668018 10.3974527
## age          0.04481205  0.1524371  0.2600621
## typeacf     -2.56775839  7.0449315 16.6576214
## typepolyp   43.84340564 51.7598345 59.6762634
```



``` r
df = full_dataset[full_dataset$type %in% c("normal", "polyp", "acf"), ]
lmm_id5 <- lme(ID5 ~ age + type, random = ~1 | patient/block/group_id, data = df,
    method = "ML")
summary(lmm_id5)
```

```
## Linear mixed-effects model fit by maximum likelihood
##   Data: df 
##        AIC      BIC   logLik
##   4821.919 4857.054 -2402.96
## 
## Random effects:
##  Formula: ~1 | patient
##         (Intercept)
## StdDev:     13.4138
## 
##  Formula: ~1 | block %in% patient
##         (Intercept)
## StdDev:    6.058864
## 
##  Formula: ~1 | group_id %in% block %in% patient
##         (Intercept) Residual
## StdDev:    10.45899  6.55551
## 
## Fixed effects:  ID5 ~ age + type 
##                 Value Std.Error  DF   t-value p-value
## (Intercept) 22.487202  6.046221 409  3.719216  0.0002
## age          0.443613  0.108064  54  4.105087  0.0001
## typeacf     -1.393703  4.575125 409 -0.304626  0.7608
## typepolyp   24.396505  4.022883 409  6.064433  0.0000
##  Correlation: 
##           (Intr) age    typecf
## age       -0.942              
## typeacf   -0.046  0.036       
## typepolyp -0.075  0.062  0.300
## 
## Standardized Within-Group Residuals:
##          Min           Q1          Med           Q3          Max 
## -1.842984214 -0.304170848 -0.004392268  0.288289489  2.829428531 
## 
## Number of Observations: 597
## Number of Groups: 
##                          patient               block %in% patient 
##                               56                              125 
## group_id %in% block %in% patient 
##                              536
```

``` r
pvalues_lst = c(pvalues_lst, summary(lmm_id5)$tTable[, "p-value"])
intervals(lmm_id5, which = "fixed")
```

```
## Approximate 95% confidence intervals
## 
##  Fixed effects:
##                   lower       est.     upper
## (Intercept)  10.6415394 22.4872019 34.332864
## age           0.2276844  0.4436132  0.659542
## typeacf     -10.3572162 -1.3937030  7.569810
## typepolyp    16.5149347 24.3965048 32.278075
```



``` r
df = full_dataset[full_dataset$type %in% c("normal", "polyp", "acf"), ]
lmm_id <- lme(indeltotal_corr ~ age + type, random = ~1 | patient/block/group_id,
    data = df, method = "ML")
summary(lmm_id)
```

```
## Linear mixed-effects model fit by maximum likelihood
##   Data: df 
##        AIC      BIC   logLik
##   6247.959 6283.094 -3115.98
## 
## Random effects:
##  Formula: ~1 | patient
##         (Intercept)
## StdDev:    30.99966
## 
##  Formula: ~1 | block %in% patient
##         (Intercept)
## StdDev:    23.68028
## 
##  Formula: ~1 | group_id %in% block %in% patient
##         (Intercept) Residual
## StdDev:    39.62731 16.42756
## 
## Fixed effects:  indeltotal_corr ~ age + type 
##                 Value Std.Error  DF   t-value p-value
## (Intercept)  69.47450  16.44729 409  4.224069  0.0000
## age           0.65271   0.29124  54  2.241146  0.0291
## typeacf      10.36663  16.56990 409  0.625630  0.5319
## typepolyp   186.22940  14.64214 409 12.718725  0.0000
##  Correlation: 
##           (Intr) age    typecf
## age       -0.945              
## typeacf   -0.073  0.058       
## typepolyp -0.120  0.100  0.288
## 
## Standardized Within-Group Residuals:
##         Min          Q1         Med          Q3         Max 
## -3.60616315 -0.14401552 -0.01902101  0.13107157  3.39058771 
## 
## Number of Observations: 597
## Number of Groups: 
##                          patient               block %in% patient 
##                               56                              125 
## group_id %in% block %in% patient 
##                              536
```

``` r
pvalues_lst = c(pvalues_lst, summary(lmm_id)$tTable[, "p-value"])
intervals(lmm_id, which = "fixed")
```

```
## Approximate 95% confidence intervals
## 
##  Fixed effects:
##                    lower       est.      upper
## (Intercept)  37.25121995  69.474496 101.697772
## age           0.07076984   0.652712   1.234654
## typeacf     -22.09686145  10.366628  42.830118
## typepolyp   157.54273469 186.229395 214.916056
```



``` r
p_corrected_polyp = p.adjust(pvalues_lst, method = "BH")
p_table_polyp = data.frame(p = pvalues_lst, p_adj = p_corrected_polyp)
```

# Whether large polyp leads to incresed sbs burden

``` r
pvalues_lst = c()
df = full_dataset[full_dataset$type %in% c("normal", "polyp", "acf", "large polyp"),
    ]
lmm_sbs18 <- lme(SBS18 ~ type + age, random = ~1 | patient/block/group_id, data = df,
    method = "ML")
summary(lmm_sbs18)
```

```
## Linear mixed-effects model fit by maximum likelihood
##   Data: df 
##        AIC     BIC    logLik
##   7877.515 7917.28 -3929.757
## 
## Random effects:
##  Formula: ~1 | patient
##         (Intercept)
## StdDev:    112.8608
## 
##  Formula: ~1 | block %in% patient
##         (Intercept)
## StdDev:    114.0591
## 
##  Formula: ~1 | group_id %in% block %in% patient
##         (Intercept) Residual
## StdDev:    105.1615 76.52254
## 
## Fixed effects:  SBS18 ~ type + age 
##                     Value Std.Error  DF   t-value p-value
## (Intercept)       -1.1590  62.67745 409 -0.018491  0.9853
## typeacf           78.0463  51.08278 409  1.527840  0.1273
## typelarge polyp 1867.8378  95.71614 409 19.514345  0.0000
## typepolyp        523.1706  48.44636 409 10.798968  0.0000
## age                5.8799   1.10574  54  5.317611  0.0000
##  Correlation: 
##                 (Intr) typecf typlrp typply
## typeacf         -0.086                     
## typelarge polyp -0.073  0.163              
## typepolyp       -0.138  0.364  0.224       
## age             -0.946  0.069  0.057  0.116
## 
## Standardized Within-Group Residuals:
##         Min          Q1         Med          Q3         Max 
## -6.68752573 -0.20795403 -0.01334227  0.20448285  3.92677618 
## 
## Number of Observations: 613
## Number of Groups: 
##                          patient               block %in% patient 
##                               56                              127 
## group_id %in% block %in% patient 
##                              539
```

``` r
pvalues_lst = c(pvalues_lst, summary(lmm_sbs18)$tTable[, "p-value"])
intervals(lmm_sbs18, which = "fixed")
```

```
## Approximate 95% confidence intervals
## 
##  Fixed effects:
##                       lower        est.       upper
## (Intercept)     -123.865624   -1.158990  121.547643
## typeacf          -21.960890   78.046305  178.053501
## typelarge polyp 1680.449714 1867.837751 2055.225789
## typepolyp        428.324913  523.170649  618.016384
## age                3.672093    5.879917    8.087741
```


``` r
df = full_dataset[full_dataset$type %in% c("normal", "polyp", "acf", "large polyp"),
    ]
lmm_sbs1 <- lme(SBS1 ~ age + type, random = ~1 | patient/block/group_id, data = df,
    method = "ML")
summary(lmm_sbs1)
```

```
## Linear mixed-effects model fit by maximum likelihood
##   Data: df 
##        AIC     BIC    logLik
##   8217.744 8257.51 -4099.872
## 
## Random effects:
##  Formula: ~1 | patient
##         (Intercept)
## StdDev:    171.4885
## 
##  Formula: ~1 | block %in% patient
##         (Intercept)
## StdDev:    163.2286
## 
##  Formula: ~1 | group_id %in% block %in% patient
##         (Intercept) Residual
## StdDev:    154.0269 80.77243
## 
## Fixed effects:  SBS1 ~ age + type 
##                     Value Std.Error  DF  t-value p-value
## (Intercept)      151.2322  91.76237 409 1.648085  0.1001
## age               12.0431   1.62113  54 7.428827  0.0000
## typeacf          106.8697  71.57543 409 1.493107  0.1362
## typelarge polyp 1286.1794 136.94142 409 9.392187  0.0000
## typepolyp        561.0451  68.32647 409 8.211241  0.0000
##  Correlation: 
##                 (Intr) age    typecf typlrp
## age             -0.945                     
## typeacf         -0.078  0.063              
## typelarge polyp -0.067  0.052  0.160       
## typepolyp       -0.126  0.106  0.363  0.222
## 
## Standardized Within-Group Residuals:
##          Min           Q1          Med           Q3          Max 
## -3.567320381 -0.242447499 -0.007418988  0.246675599  2.712814604 
## 
## Number of Observations: 613
## Number of Groups: 
##                          patient               block %in% patient 
##                               56                              127 
## group_id %in% block %in% patient 
##                              539
```

``` r
pvalues_lst = c(pvalues_lst, summary(lmm_sbs1)$tTable[, "p-value"])
intervals(lmm_sbs1, which = "fixed")
```

```
## Approximate 95% confidence intervals
## 
##  Fixed effects:
##                       lower       est.     upper
## (Intercept)      -28.415386  151.23216  330.8797
## age                8.806217   12.04311   15.2800
## typeacf          -33.256884  106.86974  246.9964
## typelarge polyp 1018.082685 1286.17942 1554.2762
## typepolyp        427.279153  561.04514  694.8111
```


``` r
df = full_dataset[full_dataset$type %in% c("normal", "polyp", "acf", "large polyp"),
    ]
lmm_sbs5 <- lme(SBS5 ~ age + type, random = ~1 | patient/block/group_id, data = df,
    method = "ML")
summary(lmm_sbs5)
```

```
## Linear mixed-effects model fit by maximum likelihood
##   Data: df 
##        AIC      BIC    logLik
##   8346.323 8386.088 -4164.161
## 
## Random effects:
##  Formula: ~1 | patient
##         (Intercept)
## StdDev:    148.2329
## 
##  Formula: ~1 | block %in% patient
##         (Intercept)
## StdDev:    173.0249
## 
##  Formula: ~1 | group_id %in% block %in% patient
##         (Intercept) Residual
## StdDev:    147.1027 119.1644
## 
## Fixed effects:  SBS5 ~ age + type 
##                     Value Std.Error  DF   t-value p-value
## (Intercept)      219.1429  88.10625 409  2.487257  0.0133
## age               16.5715   1.55011  54 10.690552  0.0000
## typeacf         -110.4743  73.49605 409 -1.503133  0.1336
## typelarge polyp 1710.7664 138.78329 409 12.326890  0.0000
## typepolyp        463.6640  70.19067 409  6.605778  0.0000
##  Correlation: 
##                 (Intr) age    typecf typlrp
## age             -0.947                     
## typeacf         -0.097  0.078              
## typelarge polyp -0.083  0.065  0.160       
## typepolyp       -0.156  0.131  0.369  0.228
## 
## Standardized Within-Group Residuals:
##         Min          Q1         Med          Q3         Max 
## -5.12231388 -0.21376856 -0.02038033  0.20923563  9.22718177 
## 
## Number of Observations: 613
## Number of Groups: 
##                          patient               block %in% patient 
##                               56                              127 
## group_id %in% block %in% patient 
##                              539
```

``` r
pvalues_lst = c(pvalues_lst, summary(lmm_sbs5)$tTable[, "p-value"])
intervals(lmm_sbs5, which = "fixed")
```

```
## Approximate 95% confidence intervals
## 
##  Fixed effects:
##                      lower       est.      upper
## (Intercept)       46.65310  219.14289  391.63268
## age               13.47645   16.57153   19.66661
## typeacf         -254.36105 -110.47434   33.41237
## typelarge polyp 1439.06375 1710.76640 1982.46906
## typepolyp        326.24835  463.66396  601.07957
```



``` r
df = full_dataset[full_dataset$type %in% c("normal", "polyp", "acf", "large polyp"),
    ]
lmm_sbs <- lme(sbstotal_corr ~ age + type, random = ~1 | patient/block/group_id,
    data = df, method = "ML")
summary(lmm_sbs)
```

```
## Linear mixed-effects model fit by maximum likelihood
##   Data: df 
##        AIC      BIC    logLik
##   9697.578 9737.344 -4839.789
## 
## Random effects:
##  Formula: ~1 | patient
##         (Intercept)
## StdDev:    156.9311
## 
##  Formula: ~1 | block %in% patient
##         (Intercept)
## StdDev:    798.5377
## 
##  Formula: ~1 | group_id %in% block %in% patient
##         (Intercept) Residual
## StdDev:    422.4994 343.0425
## 
## Fixed effects:  sbstotal_corr ~ age + type 
##                    Value Std.Error  DF   t-value p-value
## (Intercept)      381.183  265.1099 409  1.437829  0.1512
## age               36.784    4.6203  54  7.961466  0.0000
## typeacf          -24.033  219.8121 409 -0.109335  0.9130
## typelarge polyp 6728.337  463.7171 409 14.509575  0.0000
## typepolyp       1921.144  223.2051 409  8.607079  0.0000
##  Correlation: 
##                 (Intr) age    typecf typlrp
## age             -0.950                     
## typeacf         -0.169  0.140              
## typelarge polyp -0.152  0.119  0.113       
## typepolyp       -0.280  0.239  0.388  0.244
## 
## Standardized Within-Group Residuals:
##         Min          Q1         Med          Q3         Max 
## -3.81742987 -0.20692766 -0.02097916  0.19507242 10.08332883 
## 
## Number of Observations: 613
## Number of Groups: 
##                          patient               block %in% patient 
##                               56                              127 
## group_id %in% block %in% patient 
##                              539
```

``` r
pvalues_lst = c(pvalues_lst, summary(lmm_sbs)$tTable[, "p-value"])
intervals(lmm_sbs, which = "fixed")
```

```
## Approximate 95% confidence intervals
## 
##  Fixed effects:
##                      lower       est.      upper
## (Intercept)     -137.83550  381.18281  900.20111
## age               27.55907   36.78436   46.00965
## typeacf         -454.36974  -24.03316  406.30342
## typelarge polyp 5820.49634 6728.33728 7636.17823
## typepolyp       1484.16464 1921.14384 2358.12305
```

# Whether large polyp leads to incresed Indel burden

``` r
df = full_dataset[full_dataset$type %in% c("normal", "polyp", "acf", "large polyp"),
    ]
lmm_id1 <- lme(ID1 ~ age + type, random = ~1 | patient/block/group_id, data = df,
    method = "ML")
summary(lmm_id1)
```

```
## Linear mixed-effects model fit by maximum likelihood
##   Data: df 
##        AIC      BIC    logLik
##   5616.815 5656.581 -2799.408
## 
## Random effects:
##  Formula: ~1 | patient
##         (Intercept)
## StdDev:    13.73767
## 
##  Formula: ~1 | block %in% patient
##         (Intercept)
## StdDev:    17.92433
## 
##  Formula: ~1 | group_id %in% block %in% patient
##         (Intercept) Residual
## StdDev:    14.47558 14.37909
## 
## Fixed effects:  ID1 ~ age + type 
##                    Value Std.Error  DF   t-value p-value
## (Intercept)      11.9348  8.750087 409  1.363966  0.1733
## age               0.2873  0.153659  54  1.869645  0.0670
## typeacf          23.7820  7.608599 409  3.125676  0.0019
## typelarge polyp 420.2317 14.149485 409 29.699433  0.0000
## typepolyp        93.0849  7.241721 409 12.853969  0.0000
##  Correlation: 
##                 (Intr) age    typecf typlrp
## age             -0.947                     
## typeacf         -0.109  0.088              
## typelarge polyp -0.093  0.073  0.159       
## typepolyp       -0.175  0.148  0.371  0.231
## 
## Standardized Within-Group Residuals:
##         Min          Q1         Med          Q3         Max 
## -9.46860856 -0.19648049 -0.01804236  0.18295066  3.35717589 
## 
## Number of Observations: 613
## Number of Groups: 
##                          patient               block %in% patient 
##                               56                              127 
## group_id %in% block %in% patient 
##                              539
```

``` r
pvalues_lst = c(pvalues_lst, summary(lmm_id1)$tTable[, "p-value"])
intervals(lmm_id1, which = "fixed")
```

```
## Approximate 95% confidence intervals
## 
##  Fixed effects:
##                        lower        est.       upper
## (Intercept)      -5.19563922  11.9348229  29.0652850
## age              -0.01952092   0.2872876   0.5940962
## typeacf           8.88630181  23.7820180  38.6777341
## typelarge polyp 392.53055954 420.2316799 447.9328003
## typepolyp        78.90740226  93.0848651 107.2623279
```


``` r
df = full_dataset[full_dataset$type %in% c("normal", "polyp", "acf", "large polyp"),
    ]
lmm_id2 <- lme(ID2 ~ age + type, random = ~1 | patient/block/group_id, data = df,
    method = "ML")
summary(lmm_id2)
```

```
## Linear mixed-effects model fit by maximum likelihood
##   Data: df 
##        AIC      BIC    logLik
##   5324.093 5363.858 -2653.047
## 
## Random effects:
##  Formula: ~1 | patient
##         (Intercept)
## StdDev:     3.64813
## 
##  Formula: ~1 | block %in% patient
##         (Intercept)
## StdDev:    2.610739
## 
##  Formula: ~1 | group_id %in% block %in% patient
##         (Intercept) Residual
## StdDev:    14.73999 11.51002
## 
## Fixed effects:  ID2 ~ age + type 
##                     Value Std.Error  DF   t-value p-value
## (Intercept)       3.67743  3.386070 409  1.086047  0.2781
## age               0.16643  0.060615  54  2.745693  0.0082
## typeacf           6.80586  5.698243 409  1.194380  0.2330
## typelarge polyp 263.57368  9.426523 409 27.960860  0.0000
## typepolyp        53.29559  4.714019 409 11.305765  0.0000
##  Correlation: 
##                 (Intr) age    typecf typlrp
## age             -0.947                     
## typeacf         -0.129  0.099              
## typelarge polyp -0.084  0.063  0.082       
## typepolyp       -0.209  0.173  0.138  0.078
## 
## Standardized Within-Group Residuals:
##         Min          Q1         Med          Q3         Max 
## -9.87143320 -0.15665502 -0.03307176  0.12913700  3.52863212 
## 
## Number of Observations: 613
## Number of Groups: 
##                          patient               block %in% patient 
##                               56                              127 
## group_id %in% block %in% patient 
##                              539
```

``` r
pvalues_lst = c(pvalues_lst, summary(lmm_id2)$tTable[, "p-value"])
intervals(lmm_id2, which = "fixed")
```

```
## Approximate 95% confidence intervals
## 
##  Fixed effects:
##                        lower        est.       upper
## (Intercept)      -2.95163825   3.6774326  10.3065034
## age               0.04540104   0.1664302   0.2874594
## typeacf          -4.34985701   6.8058648  17.9615866
## typelarge polyp 245.11892925 263.5736818 282.0284343
## typepolyp        44.06673071  53.2955899  62.5244491
```


``` r
df = full_dataset[full_dataset$type %in% c("normal", "polyp", "acf", "large polyp"),
    ]
lmm_id5 <- lme(ID5 ~ age + type, random = ~1 | patient/block/group_id, data = df,
    method = "ML")
summary(lmm_id5)
```

```
## Linear mixed-effects model fit by maximum likelihood
##   Data: df 
##        AIC      BIC    logLik
##   4946.028 4985.793 -2464.014
## 
## Random effects:
##  Formula: ~1 | patient
##         (Intercept)
## StdDev:    13.40745
## 
##  Formula: ~1 | block %in% patient
##         (Intercept)
## StdDev:    6.087051
## 
##  Formula: ~1 | group_id %in% block %in% patient
##         (Intercept) Residual
## StdDev:    10.26773 6.870264
## 
## Fixed effects:  ID5 ~ age + type 
##                    Value Std.Error  DF   t-value p-value
## (Intercept)     22.46352  6.051653 409  3.711964  0.0002
## age              0.44398  0.108158  54  4.104919  0.0001
## typeacf         -1.47832  4.550468 409 -0.324871  0.7454
## typelarge polyp 80.55169  7.707662 409 10.450858  0.0000
## typepolyp       24.50131  3.993999 409  6.134533  0.0000
##  Correlation: 
##                 (Intr) age    typecf typlrp
## age             -0.942                     
## typeacf         -0.047  0.037              
## typelarge polyp -0.038  0.029  0.166       
## typepolyp       -0.075  0.061  0.305  0.185
## 
## Standardized Within-Group Residuals:
##          Min           Q1          Med           Q3          Max 
## -2.719187428 -0.313724377 -0.002211162  0.303094970  2.749914448 
## 
## Number of Observations: 613
## Number of Groups: 
##                          patient               block %in% patient 
##                               56                              127 
## group_id %in% block %in% patient 
##                              539
```

``` r
pvalues_lst = c(pvalues_lst, summary(lmm_id5)$tTable[, "p-value"])
intervals(lmm_id5, which = "fixed")
```

```
## Approximate 95% confidence intervals
## 
##  Fixed effects:
##                       lower       est.      upper
## (Intercept)      10.6159085 22.4635167 34.3111249
## age               0.2280227  0.4439807  0.6599387
## typeacf         -10.3869852 -1.4783166  7.4303520
## typelarge polyp  65.4620292 80.5516868 95.6413444
## typepolyp        16.6820732 24.5013144 32.3205556
```



``` r
df = full_dataset[full_dataset$type %in% c("normal", "polyp", "acf", "large polyp"),
    ]
lmm_id <- lme(indeltotal_corr ~ age + type, random = ~1 | patient/block/group_id,
    data = df, method = "ML")
summary(lmm_id)
```

```
## Linear mixed-effects model fit by maximum likelihood
##   Data: df 
##        AIC      BIC    logLik
##   6648.894 6688.659 -3315.447
## 
## Random effects:
##  Formula: ~1 | patient
##         (Intercept)
## StdDev:    28.91995
## 
##  Formula: ~1 | block %in% patient
##         (Intercept)
## StdDev:     31.7492
## 
##  Formula: ~1 | group_id %in% block %in% patient
##         (Intercept) Residual
## StdDev:    40.52617 29.61542
## 
## Fixed effects:  indeltotal_corr ~ age + type 
##                    Value Std.Error  DF   t-value p-value
## (Intercept)      70.5268  17.86130 409  3.948579  0.0001
## age               0.6333   0.31445  54  2.014046  0.0490
## typeacf           6.4598  18.61730 409  0.346981  0.7288
## typelarge polyp 793.6952  32.49262 409 24.426935  0.0000
## typepolyp       187.9040  16.75169 409 11.217017  0.0000
##  Correlation: 
##                 (Intr) age    typecf typlrp
## age             -0.947                     
## typeacf         -0.101  0.081              
## typelarge polyp -0.082  0.063  0.143       
## typepolyp       -0.164  0.138  0.315  0.180
## 
## Standardized Within-Group Residuals:
##         Min          Q1         Med          Q3         Max 
## -9.02877160 -0.18751317 -0.02829253  0.16539687  3.12950352 
## 
## Number of Observations: 613
## Number of Groups: 
##                          patient               block %in% patient 
##                               56                              127 
## group_id %in% block %in% patient 
##                              539
```

``` r
pvalues_lst = c(pvalues_lst, summary(lmm_id)$tTable[, "p-value"])
intervals(lmm_id, which = "fixed")
```

```
## Approximate 95% confidence intervals
## 
##  Fixed effects:
##                         lower        est.      upper
## (Intercept)      35.558835568  70.5267555 105.494675
## age               0.005458893   0.6333235   1.261188
## typeacf         -29.988129960   6.4598447  42.907819
## typelarge polyp 730.082877784 793.6952360 857.307594
## typepolyp       155.108445215 187.9040324 220.699620
```


``` r
p_corrected_large_polyp = p.adjust(pvalues_lst, method = "BH")
p_table_large_polyp = data.frame(p = pvalues_lst, p_adj = p_corrected_large_polyp)
```



