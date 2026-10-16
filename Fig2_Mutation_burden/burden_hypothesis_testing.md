---
title: "Burden, structural variation and telomere hypothesis testing"
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

OUT_DIR <- file.path(ROOT, "Fig2_Mutation_burden")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

FAP_CSV  <- file.path(ROOT, "data/fap_database_withexposure.csv")
CTRL_CSV <- file.path(ROOT, "data/ctrl_database_withexposure.csv")
CNV_CSV  <- file.path(ROOT, "data/CNV_SV_RT/CNV.csv")
SV_CSV   <- file.path(ROOT, "data/CNV_SV_RT/SV.csv")

# APC gene body, hg38. Used to strip the disease-driver chr5 lesion out of the
# CNV/SV counts before modelling.
APC_CHR <- "chr5"; APC_START <- 112737819; APC_END <- 112846239

RESPONSES    <- c("sbstotal_corr", "indeltotal_corr", "CNV", "SV", "RT",
                  "telomerehunter2_tel_content")
REPORT_TYPES <- c("aca", "microadenoma", "diminutive_polyp", "large_polyp")
LADDER       <- c("normal", REPORT_TYPES)

# subgroup_id separates independent APC-mutant lineages within
# one lesion, and trapped crypts from their polyp.
RANDOM_FORMULA <- ~1 | patient/block/subgroup_id
```

## Build the modelling dataset


``` r
database <- read.csv(FAP_CSV, header = TRUE, stringsAsFactors = FALSE,
                     na.strings = c("", "NA"))


group_of_sample <- setNames(database$group_id, database$sample)
DRIVER_GAP <- 1e6   # gap tolerated when reconstructing the driver segment

# Per lineage, rebuild the chr5 driver interval: seed it with the segments that
# overlap APC, then grow it over segments lying within DRIVER_GAP of it. This
# removes the initiating APC event together with the pieces it gets split into. 
apc_driver_flags <- function(seg_start, seg_end, group_id, seed) {
  drop <- rep(FALSE, length(seg_start))
  for (g in unique(group_id[!is.na(group_id)])) {
    idx <- which(group_id == g & !is.na(seg_start))
    if (!length(idx) || !any(seed[idx])) next
    s <- seg_start[idx]; e <- seg_end[idx]
    lo <- min(s[seed[idx]]); hi <- max(e[seed[idx]])
    repeat {
      near <- s <= hi + DRIVER_GAP & e >= lo - DRIVER_GAP
      nlo <- min(lo, min(s[near])); nhi <- max(hi, max(e[near]))
      if (nlo == lo && nhi == hi) break
      lo <- nlo; hi <- nhi
    }
    drop[idx[s <= hi + DRIVER_GAP & e >= lo - DRIVER_GAP]] <- TRUE
  }
  drop
}

drop_apc_cnv <- function(ev) {
  on5 <- ev$chr == APC_CHR
  seed <- on5 & ev$start <= APC_END & ev$end >= APC_START
  drop <- rep(FALSE, nrow(ev))
  drop[on5] <- apc_driver_flags(ev$start[on5], ev$end[on5], ev$group_id[on5], seed[on5])
  ev[!drop, ]
}

drop_apc_sv <- function(ev) {
  on5 <- ev$chrom1 == APC_CHR | ev$chrom2 == APC_CHR
  s1 <- ifelse(ev$chrom1 == APC_CHR, ev$pos1_start, NA_real_)
  e1 <- ifelse(ev$chrom1 == APC_CHR, ev$pos1_end,   NA_real_)
  s2 <- ifelse(ev$chrom2 == APC_CHR, ev$pos2_start, NA_real_)
  e2 <- ifelse(ev$chrom2 == APC_CHR, ev$pos2_end,   NA_real_)
  seg_s <- suppressWarnings(pmin(s1, s2, na.rm = TRUE))
  seg_e <- suppressWarnings(pmax(e1, e2, na.rm = TRUE))
  seed <- (ev$chrom1 == APC_CHR & ev$pos1_start <= APC_END & ev$pos1_end >= APC_START) |
          (ev$chrom2 == APC_CHR & ev$pos2_start <= APC_END & ev$pos2_end >= APC_START) |
          # both breakpoints on chr5 FLANKING APC: the gene lies between them
          (ev$chrom1 == APC_CHR & ev$chrom2 == APC_CHR &
           seg_s <= APC_START & seg_e >= APC_END)
  seed[is.na(seed)] <- FALSE
  drop <- rep(FALSE, nrow(ev))
  drop[on5] <- apc_driver_flags(seg_s[on5], seg_e[on5], ev$group_id[on5], seed[on5])
  ev[!drop, ]
}

cnv_events <- read.csv(CNV_CSV, header = TRUE, stringsAsFactors = FALSE)
cnv_events$group_id <- group_of_sample[cnv_events$sample]
cnv_events <- drop_apc_cnv(cnv_events)
cnv_counts <- table(cnv_events$sample)

sv_events <- read.csv(SV_CSV, header = TRUE, stringsAsFactors = FALSE)
sv_events$group_id <- group_of_sample[sv_events$sample]
sv_events <- drop_apc_sv(sv_events)
sv_counts <- table(sv_events$sample)

# Absence from the event files is a genuine zero-event call, not missingness.
m <- match(database$sample, names(cnv_counts))
database$CNV <- ifelse(is.na(m), 0L, as.integer(cnv_counts[m]))
m <- match(database$sample, names(sv_counts))
database$SV  <- ifelse(is.na(m), 0L, as.integer(sv_counts[m]))

database$type[database$type %in% c("monocryptal_adenoma", "oligocryptal_adenoma")] <- "microadenoma"
database$num_APC <- as.numeric(3 - (is.na(database$germlinehit) +
                                    is.na(database$somatic1) +
                                    is.na(database$somatic2)))
database <- database[database$coverage > 10, ]

database$type[database$num_APC < 2 &
              database$type %in% c("bc", "aca", "rc", "trapped_crypt")] <- "bystander"
drop_no_apc <- database$num_APC < 2 &
  database$type %in% c("microadenoma", "diminutive_polyp", "large_polyp")
database <- database[!drop_no_apc, ]

drop_leftover <- database$type %in% c("bc", "rc")
database <- database[!drop_leftover, ]

database$group_id[database$type == "normal"]    <- database$sample[database$type == "normal"]
database$subgroup_id[database$type == "normal"] <- database$sample[database$type == "normal"]
database$type <- relevel(factor(database$type), ref = "normal")

database$germlinehit <- "APC"

ctrl_database <- read.csv(CTRL_CSV, header = TRUE, stringsAsFactors = FALSE,
                          na.strings = c("", "NA"))
ctrl_database$num_APC <- 0
ctrl_database <- ctrl_database[ctrl_database$coverage > 10, ]
shared_cols  <- intersect(colnames(database), colnames(ctrl_database))
full_dataset <- rbind(database[, shared_cols], ctrl_database[, shared_cols])
full_dataset$germlinehit[is.na(full_dataset$germlinehit)] <- "WT"
full_dataset$germlinehit <- relevel(as.factor(full_dataset$germlinehit), ref = "WT")
full_dataset$type <- factor(full_dataset$type,
                            levels = unique(c(levels(factor(full_dataset$type)), "bystander")))
```

## Models


``` r
# The per-model summaries are verbose; capture and discard them so the knitted
# document stays readable. Assignments inside capture.output() still take
# effect here, so `results` and the fitted objects are available downstream.
invisible(capture.output({
  cat("# Whether germline APC mutation leads to increased SBS burden\n\n")
  df <- full_dataset[full_dataset$type %in% c("normal"), ]
  lmm_sbs  <- lme(sbstotal_corr ~ age:germlinehit, random = ~1 | patient/block, data = df, method = "ML")
  lmm2_sbs <- lme(sbstotal_corr ~ age,             random = ~1 | patient/block, data = df, method = "ML")
  print(anova(lmm_sbs, lmm2_sbs))
  germline_lrt_sbs <- anova(lmm_sbs, lmm2_sbs)

  cat("# SBS mutation rate in normal crypts in FAP\n\n")
  df <- database[database$type %in% c("normal"), ]
  lmm2_sbs <- lme(sbstotal_corr ~ age, random = ~1 | patient/block, data = df, method = "REML")
  print(summary(lmm2_sbs))
  print(intervals(lmm2_sbs, which = "fixed"))
  pvalues_lst <- summary(lmm2_sbs)$tTable[, "p-value"]

  cat("\n# Whether germline APC mutation leads to increased Indel burden\n\n")
  df <- full_dataset[full_dataset$type %in% c("normal"), ]
  lmm_id  <- lme(indeltotal_corr ~ age:germlinehit, random = ~1 | patient/block, data = df, method = "ML")
  lmm2_id <- lme(indeltotal_corr ~ age,             random = ~1 | patient/block, data = df, method = "ML")
  print(anova(lmm_id, lmm2_id))
  germline_lrt_id <- anova(lmm_id, lmm2_id)

  cat("# Indel mutation rate in normal crypts in FAP\n\n")
  df <- database[database$type %in% c("normal"), ]
  lmm2_id <- lme(indeltotal_corr ~ age, random = ~1 | patient/block, data = df, method = "REML")
  print(summary(lmm2_id))
  print(intervals(lmm2_id, which = "fixed"))
  pvalues_lst <- c(pvalues_lst, summary(lmm2_id)$tTable[, "p-value"])

  cat("\n# Whether germline APC mutation affects telomere content\n\n")
  df <- full_dataset[full_dataset$type %in% c("normal"), ]
  lmm_tel  <- lme(telomerehunter2_tel_content ~ age:germlinehit, random = ~1 | patient/block, data = df, method = "ML")
  lmm2_tel <- lme(telomerehunter2_tel_content ~ age,             random = ~1 | patient/block, data = df, method = "ML")
  print(anova(lmm_tel, lmm2_tel))
  germline_lrt_tel <- anova(lmm_tel, lmm2_tel)

  p_table <- data.frame(p = pvalues_lst, p_adj = p.adjust(pvalues_lst, method = "BH"))
  cat("\n## P-adjusted for SBS, ID burden and telomere content in normal crypts\n")
  print(p_table)

  ## ---------------- bystander crypts vs normal crypts ----------------
  # FAP-only (2026-09-14). This is a within-FAP contrast, but SBS, indel and
  # telomere used to be tested against a `normal` arm that pooled the 400
  # control crypts, while CNV/SV/RT already used the FAP crypts alone via a
  # FAP_ONLY switch -- two baselines inside one loop. It now matches the ladder
  # below and the signature ladder in Fig3_Signature_analysis, both FAP-only.
  # No verdict changes: SBS 0.523 -> 0.253, indel 0.283 -> 0.195,
  # telomere 0.925 -> 0.782, all non-significant before and after.
  by_fap  <- database[database$type %in% c("normal", "bystander"), ]
  by_fap$type  <- factor(as.character(by_fap$type),  levels = c("normal", "bystander"))

  by_p <- setNames(rep(NA_real_, length(RESPONSES)), RESPONSES)
  for (v in RESPONSES) {
    by_df <- by_fap
    cat("\n# Whether bystander crypts differ from normal crypts in ", v,
        "  [FAP cohort only]", "\n\n", sep = "")
    m1 <- try(lme(as.formula(paste(v, "~ age + type")),
                  random = ~1 | patient/block/subgroup_id, data = by_df, method = "ML"), silent = TRUE)
    m0 <- try(lme(as.formula(paste(v, "~ age")),
                  random = ~1 | patient/block/subgroup_id, data = by_df, method = "ML"), silent = TRUE)
    if (inherits(m1, "try-error") || inherits(m0, "try-error")) {
      cat("  MODEL FAILED; recorded as NA\n"); next
    }
    a <- anova(m1, m0)
    print(a)
    by_p[v] <- a[2, "p-value"]
  }

  # LRT p-values only; these feed the bystander rows of `results` below.
  by_results <- data.frame(response = RESPONSES,
                           p.value  = as.numeric(by_p[RESPONSES]),
                           stringsAsFactors = FALSE)

  ## ---------------- the ladder: every type against normal ----------------
  cat("crypts entering the ladder (coverage > 10; bystander is never a category):\n")
  print(table(droplevels(database$type[database$type %in% LADDER])))
  cat("\ndropped, adenomatous lesion crypts with no somatic APC mutation:",
      sum(drop_no_apc), "\n")

  rows <- list()
  for (v in RESPONSES) {
    for (k in 2:length(LADDER)) {
      category <- LADDER[k]; type_set <- LADDER[1:k]
      cat("\n\n# ", v, ", category ", category, "\n\n", sep = "")
      df <- database[database$type %in% type_set, ]
      df$type <- factor(df$type, levels = type_set)
      cat("n per type:\n"); print(table(df$type))
      forms <- list(RANDOM_FORMULA, ~1 | patient/block/group_id, ~1 | patient/block, ~1 | patient)
      m <- NULL; used <- NA_character_
      for (rf in forms) {
        m <- try(lme(as.formula(paste(v, "~ age + type")), random = rf,
                     data = df, method = "REML"), silent = TRUE)
        if (!inherits(m, "try-error")) { used <- deparse(rf); break }
        m <- NULL
      }
      cn <- paste0("type", category)
      if (is.null(m)) {
        cat("  MODEL FAILED at every random structure; recorded as NA\n")
        rows[[length(rows) + 1]] <- data.frame(
          response = v, group = category, estimate = NA_real_, std_error = NA_real_,
          ci_lower = NA_real_, ci_upper = NA_real_, p = NA_real_,
          random = NA_character_, stringsAsFactors = FALSE)
        next
      }
      if (!identical(used, deparse(RANDOM_FORMULA)))
        cat("  [random structure reduced to", used, "- fuller nesting was singular]\n")
      print(summary(m))
      ci <- try(intervals(m, which = "fixed")$fixed, silent = TRUE)
      tt <- summary(m)$tTable
      if (inherits(ci, "try-error")) { lo <- NA_real_; hi <- NA_real_ }
      else { print(ci); lo <- ci[cn, "lower"]; hi <- ci[cn, "upper"] }
      rows[[length(rows) + 1]] <- data.frame(
        response = v, group = category,
        estimate = tt[cn, "Value"], std_error = tt[cn, "Std.Error"],
        ci_lower = lo, ci_upper = hi,
        p = tt[cn, "p-value"], random = used, stringsAsFactors = FALSE)
    }
  }

  by_lrt_rows <- data.frame(
    response = by_results$response, group = "bystander",
    estimate = NA_real_, std_error = NA_real_, ci_lower = NA_real_, ci_upper = NA_real_,
    p = by_results$p.value, random = NA_character_, stringsAsFactors = FALSE)
  results <- rbind(do.call(rbind, rows), by_lrt_rows)

  ladder_rows <- is.na(results$random) == FALSE | results$group != "bystander"
  results$p_adj <- NA_real_
  results$p_adj[ladder_rows] <- p.adjust(results$p[ladder_rows], method = "BH")

  # ---------------- normal-crypt results, into the same table ----------------
  # These are the figures the Results paragraph on normal crypts quotes: the
  # age slope in FAP normal crypts and the three FAP-vs-wild-type likelihood
  # ratio tests. Until now they were printed inside capture.output() and
  # discarded, so nothing on disk held them and the manuscript's 43.5/year,
  # 21-66 and P = 0.46 could not be checked against a file. They carry their
  # OWN p_adj, from the four-value BH family in p_table above; they must not
  # join the ladder's BH family, which is why they are appended after it.
  ci_sbs <- intervals(lmm2_sbs, which = "fixed")$fixed
  ci_id  <- intervals(lmm2_id,  which = "fixed")$fixed
  tt_sbs <- summary(lmm2_sbs)$tTable
  tt_id  <- summary(lmm2_id)$tTable

  normal_rows <- rbind(
    data.frame(response = "sbstotal_corr", group = "age_slope_FAP_normal",
               estimate = tt_sbs["age", "Value"], std_error = tt_sbs["age", "Std.Error"],
               ci_lower = ci_sbs["age", "lower"], ci_upper = ci_sbs["age", "upper"],
               p = tt_sbs["age", "p-value"], random = "~1 | patient/block",
               p_adj = p_table$p_adj[2], stringsAsFactors = FALSE),
    data.frame(response = "indeltotal_corr", group = "age_slope_FAP_normal",
               estimate = tt_id["age", "Value"], std_error = tt_id["age", "Std.Error"],
               ci_lower = ci_id["age", "lower"], ci_upper = ci_id["age", "upper"],
               p = tt_id["age", "p-value"], random = "~1 | patient/block",
               p_adj = p_table$p_adj[4], stringsAsFactors = FALSE))

  lrt_row <- function(resp, a) data.frame(
    response = resp, group = "FAP_vs_wildtype_LRT", estimate = NA_real_,
    std_error = NA_real_, ci_lower = NA_real_, ci_upper = NA_real_,
    p = a[["p-value"]][2], random = "~1 | patient/block", p_adj = NA_real_,
    stringsAsFactors = FALSE)

  results <- rbind(results, normal_rows,
                   lrt_row("sbstotal_corr", germline_lrt_sbs),
                   lrt_row("indeltotal_corr", germline_lrt_id),
                   lrt_row("telomerehunter2_tel_content", germline_lrt_tel))

  cat("\n\n## BH-adjusted results, all responses\n")
}))
```

## Results table


``` r
out_csv <- file.path(OUT_DIR, "burden_and_structural_results.csv")
write.csv(results, out_csv, row.names = FALSE)

# The counts the Results paragraph on normal crypts quotes. Printed here rather
# than in the models chunk, whose output is swallowed by capture.output().
nrm <- database[database$type == "normal", ]
cat(sprintf("normal crypts: n = %d | RT events = %d | CNV events = %d | SV events = %d\n",
            nrow(nrm), sum(nrm$RT, na.rm = TRUE), sum(nrm$CNV, na.rm = TRUE),
            sum(nrm$SV, na.rm = TRUE)))
```

```
## normal crypts: n = 192 | RT events = 159 | CNV events = 10 | SV events = 10
```

``` r
knitr::kable(results, digits = 4)
```



|response                    |group                |  estimate| std_error|  ci_lower|  ci_upper|      p|random                              |  p_adj|
|:---------------------------|:--------------------|---------:|---------:|---------:|---------:|------:|:-----------------------------------|------:|
|sbstotal_corr               |aca                  |  310.1614|   60.9614|  189.8175|  430.5052| 0.0000|~1 &#124; patient/block/subgroup_id | 0.0000|
|sbstotal_corr               |microadenoma         |  721.4279|   75.4643|  572.5022|  870.3535| 0.0000|~1 &#124; patient/block/subgroup_id | 0.0000|
|sbstotal_corr               |diminutive_polyp     | 1382.0868|  127.6154| 1130.3951| 1633.7784| 0.0000|~1 &#124; patient/block/subgroup_id | 0.0000|
|sbstotal_corr               |large_polyp          | 6584.8255|  656.0106| 5290.9970| 7878.6540| 0.0000|~1 &#124; patient/block/subgroup_id | 0.0000|
|indeltotal_corr             |aca                  |   27.3723|    5.1857|   17.1352|   37.6094| 0.0000|~1 &#124; patient/block/subgroup_id | 0.0000|
|indeltotal_corr             |microadenoma         |   62.1189|    6.5455|   49.2016|   75.0362| 0.0000|~1 &#124; patient/block/subgroup_id | 0.0000|
|indeltotal_corr             |diminutive_polyp     |  152.6136|   14.0468|  124.9095|  180.3177| 0.0000|~1 &#124; patient/block/subgroup_id | 0.0000|
|indeltotal_corr             |large_polyp          |  768.0942|   41.5708|  686.1055|  850.0828| 0.0000|~1 &#124; patient/block/subgroup_id | 0.0000|
|CNV                         |aca                  |   -0.0691|    0.0445|   -0.1569|    0.0187| 0.1221|~1 &#124; patient/block             | 0.1274|
|CNV                         |microadenoma         |   -0.0618|    0.0431|   -0.1468|    0.0232| 0.1535|~1 &#124; patient/block             | 0.1535|
|CNV                         |diminutive_polyp     |    0.1447|    0.0445|    0.0568|    0.2325| 0.0014|~1 &#124; patient/block/subgroup_id | 0.0017|
|CNV                         |large_polyp          |    3.9105|    0.3535|    3.2133|    4.6078| 0.0000|~1 &#124; patient/block/subgroup_id | 0.0000|
|SV                          |aca                  |   -0.1129|    0.0419|   -0.1956|   -0.0303| 0.0076|~1 &#124; patient/block             | 0.0091|
|SV                          |microadenoma         |   -0.0625|    0.0373|   -0.1360|    0.0111| 0.0955|~1 &#124; patient/block             | 0.1042|
|SV                          |diminutive_polyp     |    0.3716|    0.0795|    0.2148|    0.5284| 0.0000|~1 &#124; patient/block/subgroup_id | 0.0000|
|SV                          |large_polyp          |    3.9548|    0.2998|    3.3634|    4.5462| 0.0000|~1 &#124; patient/block/subgroup_id | 0.0000|
|RT                          |aca                  |    0.8356|    0.4054|    0.0354|    1.6359| 0.0408|~1 &#124; patient/block/subgroup_id | 0.0466|
|RT                          |microadenoma         |    1.4106|    0.4148|    0.5921|    2.2291| 0.0008|~1 &#124; patient/block/subgroup_id | 0.0011|
|RT                          |diminutive_polyp     |    5.2428|    0.7233|    3.8163|    6.6693| 0.0000|~1 &#124; patient/block/subgroup_id | 0.0000|
|RT                          |large_polyp          |   29.5484|    2.4780|   24.6611|   34.4357| 0.0000|~1 &#124; patient/block/subgroup_id | 0.0000|
|telomerehunter2_tel_content |aca                  | -180.5791|   43.8270| -267.0980|  -94.0601| 0.0001|~1 &#124; patient/block/subgroup_id | 0.0001|
|telomerehunter2_tel_content |microadenoma         | -291.7321|   40.6104| -371.8750| -211.5892| 0.0000|~1 &#124; patient/block/subgroup_id | 0.0000|
|telomerehunter2_tel_content |diminutive_polyp     | -386.9016|   26.1450| -438.4666| -335.3367| 0.0000|~1 &#124; patient/block/subgroup_id | 0.0000|
|telomerehunter2_tel_content |large_polyp          | -450.8943|   66.6235| -582.2936| -319.4949| 0.0000|~1 &#124; patient/block/subgroup_id | 0.0000|
|sbstotal_corr               |bystander            |        NA|        NA|        NA|        NA| 0.2667|NA                                  |     NA|
|indeltotal_corr             |bystander            |        NA|        NA|        NA|        NA| 0.1156|NA                                  |     NA|
|CNV                         |bystander            |        NA|        NA|        NA|        NA| 0.3790|NA                                  |     NA|
|SV                          |bystander            |        NA|        NA|        NA|        NA| 0.3237|NA                                  |     NA|
|RT                          |bystander            |        NA|        NA|        NA|        NA| 0.0624|NA                                  |     NA|
|telomerehunter2_tel_content |bystander            |        NA|        NA|        NA|        NA| 0.6151|NA                                  |     NA|
|sbstotal_corr               |age_slope_FAP_normal |   43.5796|   10.3799|   21.1552|   66.0039| 0.0010|~1 &#124; patient/block             | 0.0042|
|indeltotal_corr             |age_slope_FAP_normal |    2.1284|    1.1342|   -0.3219|    4.5787| 0.0832|~1 &#124; patient/block             | 0.1664|
|sbstotal_corr               |FAP_vs_wildtype_LRT  |        NA|        NA|        NA|        NA| 0.4693|~1 &#124; patient/block             |     NA|
|indeltotal_corr             |FAP_vs_wildtype_LRT  |        NA|        NA|        NA|        NA| 0.1574|~1 &#124; patient/block             |     NA|
|telomerehunter2_tel_content |FAP_vs_wildtype_LRT  |        NA|        NA|        NA|        NA| 0.0789|~1 &#124; patient/block             |     NA|
