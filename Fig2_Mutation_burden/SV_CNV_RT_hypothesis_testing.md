---
title: "SV_CNV_RT"
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
pvalues_lst = c()
lmm_RT <- lme(RT ~ age + type, random = ~1 | patient/block/group_id, data = database_exclude_bystander,
    method = "ML")
summary(lmm_RT)
```

```
## Linear mixed-effects model fit by maximum likelihood
##   Data: database_exclude_bystander 
##        AIC      BIC    logLik
##   1255.082 1285.292 -618.5411
## 
## Random effects:
##  Formula: ~1 | patient
##          (Intercept)
## StdDev: 0.0001803033
## 
##  Formula: ~1 | block %in% patient
##          (Intercept)
## StdDev: 0.0003364849
## 
##  Formula: ~1 | group_id %in% block %in% patient
##         (Intercept) Residual
## StdDev:     4.78272 2.584839
## 
## Fixed effects:  RT ~ age + type 
##                    Value Std.Error  DF   t-value p-value
## (Intercept)     -1.84814 1.4355577 101 -1.287401  0.2009
## age              0.10688 0.0468572  13  2.281061  0.0400
## typeacf          0.04805 1.7563253 101  0.027360  0.9782
## typelarge polyp 36.36318 2.9673074 101 12.254606  0.0000
## typepolyp        6.62290 1.4256842 101  4.645421  0.0000
##  Correlation: 
##                 (Intr) age    typecf typlrp
## age             -0.931                     
## typeacf          0.137 -0.264              
## typelarge polyp  0.104 -0.182  0.101       
## typepolyp       -0.055 -0.086  0.133  0.081
## 
## Standardized Within-Group Residuals:
##         Min          Q1         Med          Q3         Max 
## -2.45227584 -0.18423275 -0.03466997  0.13164543  4.12453646 
## 
## Number of Observations: 212
## Number of Groups: 
##                          patient               block %in% patient 
##                               15                               34 
## group_id %in% block %in% patient 
##                              138
```

``` r
pvalues_lst = c(pvalues_lst, summary(lmm_RT)$tTable[, "p-value"])

fixed.m1 <- data.frame(fixef(lmm_RT))
intervals(lmm_RT, which = "fixed")
```

```
## Approximate 95% confidence intervals
## 
##  Fixed effects:
##                        lower        est.      upper
## (Intercept)     -4.662115843 -1.84813793  0.9658400
## age              0.006856142  0.10688402  0.2069119
## typeacf         -3.394693379  0.04805265  3.4907987
## typelarge polyp 30.546672382 36.36318343 42.1796945
## typepolyp        3.828279232  6.62290305  9.4175269
```


``` r
lmm_SV <- lme(SV ~ age + type, random = ~1 | patient/block/group_id, data = database_exclude_bystander,
    method = "ML")
summary(lmm_SV)
```

```
## Linear mixed-effects model fit by maximum likelihood
##   Data: database_exclude_bystander 
##        AIC      BIC    logLik
##   488.0698 518.2791 -235.0349
## 
## Random effects:
##  Formula: ~1 | patient
##          (Intercept)
## StdDev: 5.287292e-05
## 
##  Formula: ~1 | block %in% patient
##         (Intercept)
## StdDev:   0.9479009
## 
##  Formula: ~1 | group_id %in% block %in% patient
##         (Intercept)  Residual
## StdDev:   0.2577352 0.5539347
## 
## Fixed effects:  SV ~ age + type 
##                     Value Std.Error  DF   t-value p-value
## (Intercept)     -0.942189 0.5833125 101 -1.615240  0.1094
## age              0.033850 0.0187862  13  1.801884  0.0948
## typeacf         -0.356371 0.2173127 101 -1.639899  0.1041
## typelarge polyp  4.316418 0.4621791 101  9.339276  0.0000
## typepolyp        0.738746 0.2265543 101  3.260789  0.0015
##  Correlation: 
##                 (Intr) age    typecf typlrp
## age             -0.945                     
## typeacf          0.004 -0.097              
## typelarge polyp  0.108 -0.212  0.201       
## typepolyp       -0.104 -0.019  0.472  0.352
## 
## Standardized Within-Group Residuals:
##         Min          Q1         Med          Q3         Max 
## -5.01366537 -0.19420390 -0.02537418  0.06619652  4.01266947 
## 
## Number of Observations: 212
## Number of Groups: 
##                          patient               block %in% patient 
##                               15                               34 
## group_id %in% block %in% patient 
##                              138
```

``` r
pvalues_lst = c(pvalues_lst, summary(lmm_SV)$tTable[, "p-value"])

fixed.m1 <- data.frame(fixef(lmm_SV))
intervals(lmm_SV, which = "fixed")
```

```
## Approximate 95% confidence intervals
## 
##  Fixed effects:
##                        lower        est.      upper
## (Intercept)     -2.085597789 -0.94218949 0.20121881
## age             -0.006253104  0.03385049 0.07395408
## typeacf         -0.782346958 -0.35637094 0.06960509
## typelarge polyp  3.410454982  4.31641763 5.22238028
## typepolyp        0.294654438  0.73874582 1.18283719
```


``` r
lmm_CNV <- lme(CNV ~ age + type, random = ~1 | patient/block/group_id, data = database_exclude_bystander,
    method = "ML")
summary(lmm_CNV)
```

```
## Linear mixed-effects model fit by maximum likelihood
##   Data: database_exclude_bystander 
##        AIC      BIC    logLik
##   445.2849 475.4942 -213.6425
## 
## Random effects:
##  Formula: ~1 | patient
##          (Intercept)
## StdDev: 9.060635e-06
## 
##  Formula: ~1 | block %in% patient
##         (Intercept)
## StdDev:     2.05331
## 
##  Formula: ~1 | group_id %in% block %in% patient
##         (Intercept)  Residual
## StdDev: 1.66069e-05 0.4586464
## 
## Fixed effects:  CNV ~ age + type 
##                     Value Std.Error  DF   t-value p-value
## (Intercept)     -1.463949 1.1873453 101 -1.232960  0.2205
## age              0.058181 0.0379924  13  1.531380  0.1496
## typeacf         -0.001952 0.1566347 101 -0.012464  0.9901
## typelarge polyp  5.005151 0.3487107 101 14.353305  0.0000
## typepolyp        0.562117 0.1742789 101  3.225388  0.0017
##  Correlation: 
##                 (Intr) age    typecf typlrp
## age             -0.952                     
## typeacf         -0.003 -0.039              
## typelarge polyp  0.033 -0.078  0.284       
## typepolyp       -0.032 -0.020  0.579  0.484
## 
## Standardized Within-Group Residuals:
##          Min           Q1          Med           Q3          Max 
## -5.857985134 -0.033348889 -0.005186997  0.009295735  5.043658279 
## 
## Number of Observations: 212
## Number of Groups: 
##                          patient               block %in% patient 
##                               15                               34 
## group_id %in% block %in% patient 
##                              138
```

``` r
pvalues_lst = c(pvalues_lst, summary(lmm_CNV)$tTable[, "p-value"])

fixed.m1 <- data.frame(fixef(lmm_CNV))
intervals(lmm_CNV, which = "fixed")
```

```
## Approximate 95% confidence intervals
## 
##  Fixed effects:
##                       lower         est.     upper
## (Intercept)     -3.79138136 -1.463948972 0.8634834
## age             -0.02292309  0.058180804 0.1392847
## typeacf         -0.30898739 -0.001952217 0.3050830
## typelarge polyp  4.32160896  5.005151137 5.6886933
## typepolyp        0.22049586  0.562117073 0.9037383
```


``` r
p_corrected = p.adjust(pvalues_lst, method = "BH")
p_table = data.frame(p = pvalues_lst, p_adj = p_corrected)
```
