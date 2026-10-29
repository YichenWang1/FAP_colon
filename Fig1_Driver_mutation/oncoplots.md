---
title: "Driver mutation oncoplot (combined)"
author: "Yichen Wang"
output:
  html_document:
    keep_md: yes
---




## Setup


``` r
suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
  library(RColorBrewer)
})

NO_BORDER <- theme(plot.background = element_rect(fill = "white", colour = NA))
theme_set(theme_get() + NO_BORDER)
```


``` r
ROOT     <- "/lustre/scratch125/casm/teams/team267/users/yw2/colon/github_repo"
OUT_DIR  <- file.path(ROOT, "Fig1_Driver_mutation")
FAP_CSV   <- file.path(ROOT, "data/fap_database_withexposure.csv")
CNV_CSV  <- file.path(ROOT, "data/CNV_SV_RT/CNV.csv")

CATEGORIES <- c("bc", "aca", "microadenoma", "diminutive", "large", "trapped")
SEVERITY_ORDER <- c("bc", "aca", "microadenoma", "diminutive", "large")
CATEGORY_TITLES <- c(
  bc      = "Bifurcating crypts (non-trapped)",
  aca     = "Aberrant Crypt Aggregates",
  microadenoma = "Microadenomas",
  diminutive = "Diminutive polyps",
  large   = "Large polyps",
  trapped = "Trapped bystanders (normal/RC)"
)

VC_LEVELS <- c("Frame_Shift_Del", "Missense_Mutation", "Nonsense_Mutation", "Multi_Hit",
               "Frame_Shift_Ins", "Splice_Site", "LOH", "None")
vc_cols <- c(RColorBrewer::brewer.pal(n = length(VC_LEVELS) - 1, name = "Paired"), "#EEEEEE")
names(vc_cols) <- VC_LEVELS

CLASS_MAP <- c(Nonsense = "Nonsense_Mutation", Missense = "Missense_Mutation",
               Splice = "Splice_Site", Frame_Shift_Del = "Frame_Shift_Del", Frame_Shift_Ins = "Frame_Shift_Ins")
```

## Category membership and driver calls


``` r
## ---------------- category membership ----------------
database <- read.csv(FAP_CSV, header = TRUE, stringsAsFactors = FALSE, check.names = FALSE, na.strings = c("", "NA"))
database$is_trapped <- !is.na(database$trapped_in_polyp)
database$is_cnv <- database$somatic1 %in% c("Chr5LOH", "Chr5CNNLOH") | database$somatic2 %in% c("Chr5LOH", "Chr5CNNLOH")

category_rows <- list(
  bc      = database[database$type == "bc", ],
  aca     = database[database$type == "aca", ],
  microadenoma = database[database$type %in% c("monocryptal_adenoma", "oligocryptal_adenoma"), ],
  diminutive = database[database$type == "diminutive_polyp", ],
  large   = database[database$type == "large_polyp", ],
  trapped = database[database$type == "trapped_crypt", ]
)
```


``` r
somatic_drivers <- read.csv(file.path(ROOT, "data/driver_mutation/somatic_drivers.csv"), stringsAsFactors = FALSE)
somatic_drivers <- somatic_drivers[somatic_drivers$is.supressor == "Y" | somatic_drivers$is.hotspot == "Y", ]
trapped_by_group <- tapply(database$is_trapped, database$group_id, any) 
somatic_drivers$is_trapped <- trapped_by_group[somatic_drivers$Group_ID]

## Second hits: a non-APC driver mutation that falls inside a copy number
## LOSS or LOH segment (minor == 0) carried by at least one of the same crypts
## gets an LOH or Loss event (the segment's CNV.status in the shared crypts)
## for that gene, as chr5 LOH does for APC (build_events).
cnv_all <- read.csv(CNV_CSV, stringsAsFactors = FALSE)

second_hits <- local({
  sd2 <- somatic_drivers[somatic_drivers$Gene != "APC", ]
  sd2$chr <- ifelse(grepl("^chr", sd2$Chr), sd2$Chr, paste0("chr", sd2$Chr))
  hits <- data.frame(Gene = character(0), Group_ID = character(0), Variant_Classification = character(0),
                     stringsAsFactors = FALSE)
  for (i in seq_len(nrow(sd2))) {
    r <- sd2[i, ]
    o <- cnv_all[cnv_all$group_id == r$Group_ID & cnv_all$chr == r$chr &
                 cnv_all$minor == 0 &
                 cnv_all$start <= as.numeric(r$Start) &
                 cnv_all$end   >= as.numeric(r$Start), ]
    if (!nrow(o)) next
    shared <- intersect(trimws(strsplit(r$Sample, ";")[[1]]), unique(o$sample))
    if (!length(shared)) next
    cat(sprintf("second hit: %-8s %-13s in %-18s (%s, %d/%d crypts shared)\n",
                r$Gene, r$Mutation, r$Group_ID,
                paste(unique(o$CNV.status), collapse = ","), length(shared),
                length(trimws(strsplit(r$Sample, ";")[[1]]))))
    hits <- rbind(hits, data.frame(Gene = r$Gene, Group_ID = r$Group_ID,
                                   Variant_Classification = unique(o$CNV.status[o$sample %in% shared]),
                                   stringsAsFactors = FALSE))
  }
  unique(hits)
})
```

```
## second hit: TP53     R196*         in PD42778J_PLP_09    (Loss, 3/3 crypts shared)
## second hit: SOX9     Y319*         in PD44721I_PLP_02_1  (LOH, 1/4 crypts shared)
## second hit: SMAD4    R361H         in PD42778H_PLP_01    (Loss, 10/10 crypts shared)
```

``` r
driver_details <- list(
  bc      = somatic_drivers[somatic_drivers$Type == "bc" & !somatic_drivers$is_trapped, ],
  aca     = somatic_drivers[somatic_drivers$Type == "aca" & !somatic_drivers$is_trapped, ],
  microadenoma = somatic_drivers[somatic_drivers$Type %in% c("monocryptal_adenoma", "oligocryptal_adenoma"), ],
  diminutive = somatic_drivers[somatic_drivers$Type == "diminutive_polyp", ],
  large   = somatic_drivers[somatic_drivers$Type == "large_polyp", ],
  trapped = somatic_drivers[somatic_drivers$Type == "trapped_crypt", ]
)
```


``` r
## ---------------- build (Gene, Group_ID, Variant_Classification) long table per category ----------------
build_events <- function(cat) {
  muts <- driver_details[[cat]]
  ev <- data.frame(Gene = muts$Gene, Group_ID = muts$Group_ID,
                    Variant_Classification = unname(CLASS_MAP[muts$Mutation_Type]), stringsAsFactors = FALSE)
  # copy number events: chr5 LOH for APC, and overlapping LOH/Loss for other
  # drivers mutated in this category (second_hits above)
  cnv_groups <- unique(category_rows[[cat]]$group_id[category_rows[[cat]]$is_cnv])
  ev <- rbind(ev,
              data.frame(Gene = rep("APC", length(cnv_groups)), Group_ID = cnv_groups,
                         Variant_Classification = rep("LOH", length(cnv_groups)), stringsAsFactors = FALSE),
              second_hits[paste(second_hits$Gene, second_hits$Group_ID) %in% paste(ev$Gene, ev$Group_ID), ])
  # collapse >1 distinct classification to Multi_Hit
  aggregate(Variant_Classification ~ Gene + Group_ID, data = ev,
            FUN = function(x) if (length(unique(x)) > 1) "Multi_Hit" else x[1])
}

## order genes by frequency (desc) so the most-mutated gene's samples cluster leftmost.
memo_sort_samples <- function(bin_mat) {
  gene_order <- order(rowSums(bin_mat), decreasing = TRUE)
  m <- bin_mat[gene_order, , drop = FALSE]
  weights <- 2^(rev(seq_len(nrow(m))) - 1)
  score <- colSums(m * weights)
  list(gene_order = rownames(m), sample_order = colnames(m)[order(-score)])
}
```

## Combined oncoplot (Aberrant Crypt Aggregates | Microadenomas | Diminutive polyps | Large polyps).


``` r
COMBINED_ORDER <- c("aca", "microadenoma", "diminutive", "large")

all_groups_by_cat <- lapply(COMBINED_ORDER, function(cat) unique(category_rows[[cat]]$group_id))
names(all_groups_by_cat) <- COMBINED_ORDER

ev_all <- do.call(rbind, lapply(COMBINED_ORDER, build_events))
genes_all <- unique(ev_all$Gene)
full_groups_all <- unlist(all_groups_by_cat, use.names = FALSE)
n_total_all <- length(full_groups_all)

bin_mat_all <- matrix(0, nrow = length(genes_all), ncol = n_total_all,
                       dimnames = list(genes_all, full_groups_all))
for (i in seq_len(nrow(ev_all))) bin_mat_all[ev_all$Gene[i], ev_all$Group_ID[i]] <- 1

# order: block by category (aca, microadenoma, diminutive, large, so the
# adenoma subtypes stay contiguous and in that order), mutated
# samples within each block, WT samples appended per block.
sample_order_all <- unlist(lapply(COMBINED_ORDER, function(cat) {
  grp <- all_groups_by_cat[[cat]]
  mutated <- grp[colSums(bin_mat_all[, grp, drop = FALSE]) > 0]
  wt <- setdiff(grp, mutated)
  if (length(mutated) > 0) {
    ord <- memo_sort_samples(bin_mat_all[, mutated, drop = FALSE])
    c(ord$sample_order, sort(wt))
  } else sort(wt)
}))
# Genes are ordered by the number of lesions carrying a mutation, which is the
# same quantity as the per-panel percentages printed beside each row. Ties are
# broken by how many of those lesions are Multi_Hit
n_lesions_all <- rowSums(bin_mat_all)
n_multihit_all <- tapply(ev_all$Variant_Classification == "Multi_Hit", ev_all$Gene, sum)
n_multihit_all <- as.integer(n_multihit_all[names(n_lesions_all)])
n_multihit_all[is.na(n_multihit_all)] <- 0L
gene_order_all <- names(n_lesions_all)[order(-n_lesions_all, -n_multihit_all)]

plot_df <- expand.grid(Gene = gene_order_all, Group_ID = sample_order_all, stringsAsFactors = FALSE)
plot_df <- merge(plot_df, ev_all, by = c("Gene", "Group_ID"), all.x = TRUE)
plot_df$Variant_Classification[is.na(plot_df$Variant_Classification)] <- "None"
plot_df$Variant_Classification <- factor(plot_df$Variant_Classification, levels = VC_LEVELS)
plot_df$Gene <- factor(plot_df$Gene, levels = rev(gene_order_all))
plot_df$Group_ID <- factor(plot_df$Group_ID, levels = sample_order_all)

pct_by_type <- do.call(rbind, lapply(COMBINED_ORDER, function(cat) {
  grp <- all_groups_by_cat[[cat]]
  data.frame(Gene = factor(gene_order_all, levels = rev(gene_order_all)),
             cat = cat,
             pct = sapply(gene_order_all, function(g) round(sum(bin_mat_all[g, grp] > 0) / length(grp) * 100)),
             stringsAsFactors = FALSE)
}))

# Four visual panels, each with its own strip title/N.
cat_lookup <- setNames(rep(COMBINED_ORDER, sapply(all_groups_by_cat, length)), full_groups_all)

panel_labels <- c(
  aca          = paste0(CATEGORY_TITLES[["aca"]], "\n(N = ", length(all_groups_by_cat$aca), ")"),
  microadenoma = paste0("Microadenomas\n(N = ", length(all_groups_by_cat$microadenoma), ")"),
  diminutive   = paste0("Diminutive polyps\n(N = ", length(all_groups_by_cat$diminutive), ")"),
  large        = paste0("Large polyps\n(N = ", length(all_groups_by_cat$large), ")")
)
plot_df$Category <- factor(panel_labels[cat_lookup[as.character(plot_df$Group_ID)]],
                            levels = panel_labels[COMBINED_ORDER])
plot_df$Gene <- factor(as.character(plot_df$Gene), levels = rev(gene_order_all))
pct_by_type$Gene <- factor(as.character(pct_by_type$Gene), levels = rev(gene_order_all))
pct_by_type$Category <- factor(panel_labels[pct_by_type$cat], levels = panel_labels[COMBINED_ORDER])

pct_max <- max(pct_by_type$pct, na.rm = TRUE)
n_cat <- lengths(all_groups_by_cat)[COMBINED_ORDER]

FONT_FAMILY <- "Helvetica"
fig_width   <- 7.09   # ~180mm, full two-column journal figure width
fig_height  <- 2.70   # ~69mm
legend_height <- 0.22
```



``` r
blocks <- vector("list", 2 * length(COMBINED_ORDER))
widths <- numeric(2 * length(COMBINED_ORDER))
for (i in seq_along(COMBINED_ORDER)) {
  cat_i <- COMBINED_ORDER[i]
  grp <- all_groups_by_cat[[cat_i]]
  df_cat <- plot_df[as.character(plot_df$Group_ID) %in% grp, ]
  pct_cat <- pct_by_type[pct_by_type$cat == cat_i, ]

  tile_cat <- ggplot(df_cat, aes(x = Group_ID, y = Gene, fill = Variant_Classification)) +
    geom_tile(color = "white", linewidth = 0.15) +
    scale_fill_manual(values = vc_cols, breaks = setdiff(VC_LEVELS, "None"), name = NULL, drop = FALSE) +
    labs(title = panel_labels[[cat_i]]) +
    theme_minimal(base_family = FONT_FAMILY) +
    theme(axis.text.x = element_blank(), axis.ticks = element_blank(), axis.title = element_blank(),
          axis.text.y = if (i == 1) element_text(size = 5.5, colour = "black", face = "italic") else element_blank(),
          panel.grid = element_blank(), plot.margin = margin(1, 1, 1, 1),
          legend.position = "none",
          plot.title = element_text(face = "bold", size = 7, hjust = 0.5))

  pct_plot_cat <- ggplot(pct_cat, aes(x = pct, y = Gene)) +
    geom_col(fill = "grey40", width = 0.7, na.rm = TRUE) +
    geom_text(aes(label = paste0(pct, "%")), hjust = -0.15, size = 5.5 / .pt,
              family = FONT_FAMILY, na.rm = TRUE) +
    scale_x_continuous(limits = c(0, pct_max * 2.1), expand = c(0, 0)) +
    scale_y_discrete(drop = FALSE) +
    coord_cartesian(clip = "off") +
    theme_void() +
    theme(plot.margin = margin(1, if (i == length(COMBINED_ORDER)) 10 else 6, 1, 1))

  blocks[[2 * i - 1]] <- tile_cat
  blocks[[2 * i]]     <- pct_plot_cat
  # tile width scales with that category's own sample count (so column width
  # stays visually consistent across panels)
  widths[2 * i - 1] <- n_cat[[cat_i]]
  widths[2 * i]     <- 6
}

main_row <- patchwork::wrap_plots(blocks, nrow = 1, widths = widths)

legend_ref <- ggplot(plot_df, aes(x = Group_ID, y = Gene, fill = Variant_Classification)) +
  geom_tile() +
  scale_fill_manual(values = vc_cols, breaks = setdiff(VC_LEVELS, "None"), name = NULL, drop = FALSE) +
  guides(fill = guide_legend(nrow = 1)) +
  theme(text = element_text(family = FONT_FAMILY),
        legend.position = "bottom",
        legend.text = element_text(size = 5.5, margin = margin(l = 1, unit = "pt")),
        legend.key.spacing.x = unit(3, "mm"),
        legend.key.size = unit(2.5, "mm"))
legend_grob <- ggplotGrob(legend_ref)$grobs[[which(sapply(ggplotGrob(legend_ref)$grobs, `[[`, "name") == "guide-box")]]

layout <- main_row / patchwork::wrap_elements(full = legend_grob) +
  patchwork::plot_layout(heights = c(fig_height - legend_height, legend_height))

out_pdf <- file.path(OUT_DIR, "combined_oncoplot.pdf")
ggsave(out_pdf, layout, width = fig_width, height = fig_height, limitsize = FALSE)

layout
```

![](oncoplots_files/figure-html/combined-plot-1.png)<!-- -->
