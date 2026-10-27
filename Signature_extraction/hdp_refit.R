library(RColorBrewer)
library(ggtree)
library(lsa)
library(lattice)
library(sigfit)
library(ape)

# --- Path configuration ---
HDP_DIR       <- "/lustre/scratch125/casm/teams/team267/users/yw2/colon/HDP/snv"
SBS_INPUT     <- "/lustre/scratch125/casm/teams/team267/users/yw2/colon/LCMB_analysis/sbs_fap_control_perbranch.hdp.txt"
FAP_SNP_DIR   <- "/lustre/scratch125/casm/teams/team267/users/yw2/colon/LCMB_analysis/FAP/snp"
NORM_SNP_DIR  <- "/lustre/scratch125/casm/teams/team267/users/yw2/colon/LCMB_analysis/normal/snp"
DATABASE      <- "/lustre/scratch125/casm/teams/team267/users/yw2/colon/metadata/fap_database_withexposure.csv"
HENRY_SIGS    <- "/lustre/scratch125/casm/teams/team267/users/yw2/public/Henry/SBS_novel_sig_trinucs.txt"
FAP_PATIENTS  <- "/lustre/scratch125/casm/teams/team267/users/yw2/colon/LCMB_analysis/FAP/patients.txt"
NORM_PATIENTS <- "/lustre/scratch125/casm/teams/team267/users/yw2/colon/LCMB_analysis/normal/patients.txt"

dir.create(file.path(HDP_DIR, "refit"),                showWarnings=FALSE, recursive=TRUE)
dir.create(file.path(HDP_DIR, "plots"),                showWarnings=FALSE, recursive=TRUE)
dir.create(file.path(HDP_DIR, "tree_plots", "fap"),    showWarnings=FALSE, recursive=TRUE)
dir.create(file.path(HDP_DIR, "tree_plots", "control"),showWarnings=FALSE, recursive=TRUE)

# Load reference signatures
data("cosmic_signatures_v3.2")
features <- paste0(substr(colnames(cosmic_signatures_v3.2),1,1),'[',substr(colnames(cosmic_signatures_v3.2),2,2),'>',substr(colnames(cosmic_signatures_v3.2),6,6),']',substr(colnames(cosmic_signatures_v3.2),3,3))
ref<-t(cosmic_signatures_v3.2)
rownames(ref)<-features

mut.cols = rep(c("dodgerblue","black","red","grey70","olivedrab3","plum2"),each=16)

# Load HDP signatures (CSV: rows = trinucleotides, cols = X0-X9)
hdp_sigs <- read.csv(file.path(HDP_DIR, "combined_HDP_sigs.csv"), row.names=1, check.names=F)
rownames(hdp_sigs) <- features
colnames(hdp_sigs) <- paste0('Signature', 0:(ncol(hdp_sigs)-1))
rownames(hdp_sigs) <- rownames(ref)
hdp_sigs <- hdp_sigs[,-1]   # Remove Signature0 (noise component)

# Load HDP exposures (CSV: rows = samples, cols = 0-9)
hdp_exposures <- read.csv(file.path(HDP_DIR, "combined_HDP_exposure.csv"), row.names=1, check.names=F)
hdp_exposures <- hdp_exposures[,-1]   # Remove Component 0
colnames(hdp_exposures) <- paste0('Signature', 1:ncol(hdp_exposures))

input_for_hdp <- read.table(SBS_INPUT, check.names=F)
input_for_hdp <- input_for_hdp[apply(input_for_hdp,1,sum)>50,]

hdp_exposures <- hdp_exposures[rownames(input_for_hdp), , drop=FALSE]


# ---------Assess cosine similarities for all reference signatures-----------
cosine_matrix=data.frame(matrix(nrow=ncol(hdp_sigs), ncol=ncol(ref)))
rownames(cosine_matrix)=colnames(hdp_sigs)
colnames(cosine_matrix)=colnames(ref)

for (n in 1:nrow(cosine_matrix)) {
  for (m in 1:ncol(cosine_matrix)) {
    cosine_matrix[n,m] <- cosine(x=hdp_sigs[,rownames(cosine_matrix)[n]],
                                 y=ref[,colnames(cosine_matrix)[m]])
  }
}

write.table(cosine_matrix, file.path(HDP_DIR, "Cosine_similarities_cosmic.txt"), sep="\t", quote=F)

# Plot output
color.palette = colorRampPalette(c("white", "orange", "purple"))
pdf(file.path(HDP_DIR, "cosine_similarity.pdf"), width=9, height=4)
levelplot(t(cosine_matrix[dim(cosine_matrix)[1]:1,]),col.regions=color.palette, aspect="fill", scales=list(x=list(rot=90)),
          colorkey=list(title="Cosine similarity", title.gpar=list(cex=0.9)))
dev.off()

# ----------Assess cosine similarities for all suspected sigs -------------
gdsigs=c('SBS1', 'SBS2','SBS5', 'SBS13', 'SBS17a', 'SBS17b', 'SBS18', 'SBS28', 'SBS31', 'SBS35', 'SBS88', 'SBS89','SBS93')
# gdsigs=c("SBS1","SBS2","SBS5","SBS13","SBS17b","SBS18","SBS35","SBS40","SBS88","SBS89")
Henry_sigs=read.table(HENRY_SIGS, header=T)

ref<-data.frame(t(cosmic_signatures_v3.2))
rownames(ref)<-features
ref=cbind(ref[,gdsigs],Henry_sigs$SBSC)
colnames(ref)[ncol(ref)] = 'SBS112'   # Henry's SBSC is our SBS112
gdsigs=c(gdsigs,"SBS112")
ref=ref[,gdsigs]

cosine_matrix=data.frame(matrix(nrow=ncol(hdp_sigs), ncol=ncol(ref)))
rownames(cosine_matrix)=colnames(hdp_sigs)
colnames(cosine_matrix)=colnames(ref)

for (n in 1:nrow(cosine_matrix)) {
  for (m in 1:ncol(cosine_matrix)) {
    cosine_matrix[n,m] <- cosine(x=hdp_sigs[,rownames(cosine_matrix)[n]],
                                 y=ref[,colnames(cosine_matrix)[m]])
  }
}

write.table(cosine_matrix, file.path(HDP_DIR, "Cosine_similarities_ref.txt"), sep="\t", quote=F)


# Plot output. levelplot as before, with three changes: APOBEC dropped (it is
# (SBS2 + SBS13)/2, built above only for this comparison, not a COSMIC
# signature), Nature two-column sizing, and the colour-key title placed ABOVE
# the key via colorkey$title.control$side (lattice >= 0.21) rather than rotated
# alongside it. Matches script/extended_figure_plotting/replot_cosine_similarity_ref.R.
color.palette = colorRampPalette(c("white", "orange", "purple"))
cos_m <- cosine_matrix[, setdiff(colnames(cosine_matrix), "APOBEC"), drop=FALSE]
pdf(file.path(HDP_DIR, "cosine_similarity_ref.pdf"), width=180/25.4, height=48/25.4)
print(levelplot(t(cos_m[dim(cos_m)[1]:1,]), col.regions=color.palette, aspect="fill",
                xlab=NULL, ylab=NULL,
                scales=list(x=list(rot=90), tck=c(0,0)),
                colorkey=list(space="right", width=0.9,
                              title="Cosine\nsimilarity",
                              title.control=list(side="top", padding=1)),
                par.settings=list(fontsize=list(text=5.8),
                                  layout.widths=list(right.padding=6))))
dev.off()


ref[,'SBS112']=hdp_sigs$Signature8

# --------------First iteration; decomposed hdp sigs into all suspected signatures---------
sigs_to_decompose=rowSums(cosine_matrix>=0.9)
for(n in 1:nrow(cosine_matrix)){
  print(paste0(rownames(cosine_matrix)[n],": ",paste(colnames(cosine_matrix)[cosine_matrix[n,]>=0.9],collapse=",")))
}
# Select which signatures to decompose into reference signatures
sigs_to_deconv=names(sigs_to_decompose)[sigs_to_decompose!=1]
sigs_to_deconv

# Hard-code the dirty SBS1 signatures(HDP Signature 1) to be decomposed
sigs_to_deconv=c(sigs_to_deconv,'Signature1')


signatures=t(ref[,gdsigs])
signatures = signatures[rownames(signatures)!='SBS112',]
signature_fractionR1 = matrix(NA,nrow=nrow(signatures),ncol=length(sigs_to_deconv))
rownames(signature_fractionR1) = rownames(signatures)
colnames(signature_fractionR1) = sigs_to_deconv
maxiter <- 1000

for (j in 1:length(sigs_to_deconv)) {
  freqs = hdp_sigs[,sigs_to_deconv[j]]
  freqs[is.na(freqs)] = 0
  # EM algorithm to estimate the signature contribution
  alpha = runif(nrow(signatures)); alpha=alpha/sum(alpha) # Random start (seems to give ~identical results)
  for (iter in 1:maxiter) {
    contr = t(array(alpha,dim=c(nrow(signatures),96))) * t(signatures)
    probs = contr/array(rowSums(contr),dim=dim(contr))
    probs = probs * freqs
    probs[is.na(probs)] = 0
    old_alpha = alpha
    alpha = colSums(probs)/sum(probs)
    if (sum(abs(alpha-old_alpha))<1e-5) {
      break
    }
  }
  # Saving the signature contributions for the sample
  signature_fractionR1[,j] = alpha
}

# Plot initial deconvolution
pdf(file.path(HDP_DIR, "sigs_decomposition_r1.pdf"), width=9, height=4)
color.palette = colorRampPalette(c("white", "orange", "purple"))
levelplot((signature_fractionR1[nrow(signature_fractionR1):1,]),col.regions=color.palette, aspect="fill", scales=list(x=list(rot=90)))
dev.off()


sigs_deconv_R2=list()
for(n in colnames(hdp_sigs)){
  if (n %in% sigs_to_deconv){
    sigs_deconv_R2[[n]]=rownames(signature_fractionR1)[signature_fractionR1[,n]>0.1]
  }
  else{
    sigs_deconv_R2[[n]]=colnames(cosine_matrix)[cosine_matrix[n,]==max(cosine_matrix[n,])]
  }

}

sigs_deconv_R2

# Combine decomposed and undecomposed signatures
ref_sigs_R2=sort(unique(unlist(sigs_deconv_R2)))
signature_fractionR2=matrix(0,ncol=length(colnames(hdp_sigs)),nrow=length(ref_sigs_R2))
rownames(signature_fractionR2)=ref_sigs_R2
colnames(signature_fractionR2)=colnames(hdp_sigs)

pdf(file.path(HDP_DIR, "HDPsig_decomposed.pdf"), width=9, height=4)
# Repeat the deconvolution with the identified constitutive signatures
for(s in colnames(hdp_sigs)){
  gdsigs <- sigs_deconv_R2[[s]]
  signatures <- t(ref[,gdsigs])
  signatures = signatures[rownames(signatures)!='SBS112',]
  if (length(gdsigs)==1){
    signature_fractionR2[gdsigs,s] = 1
    next
  }
  signature_fraction = matrix(NA,nrow=nrow(signatures),ncol=length(colnames(hdp_sigs)))
  rownames(signature_fraction) = rownames(signatures)
  colnames(signature_fraction) = colnames(hdp_sigs)
  maxiter <- 1000

  freqs = hdp_sigs[,s]
  freqs[is.na(freqs)] = 0

  alpha = runif(nrow(signatures)); alpha=alpha/sum(alpha) # Random start (seems to give ~identical results)
  for (iter in 1:maxiter) {
    contr = t(array(alpha,dim=c(nrow(signatures),96))) * t(signatures)
    probs = contr/array(rowSums(contr),dim=dim(contr))
    probs = probs * freqs
    probs[is.na(probs)] = 0
    old_alpha = alpha
    alpha = colSums(probs)/sum(probs)
    if (sum(abs(alpha-old_alpha))<1e-5) {
      break
    }
  }
  # Saving the signature contributions for the sample
  signature_fractionR2[gdsigs,s] = alpha
  reconsbs <- rep(0,96)
  for (g in gdsigs) {
    reconsbs=reconsbs+(ref[,g]*alpha[g])
  }
  cosine_reconst=cosine(x=reconsbs, y=hdp_sigs[,s])
  print(paste0(s,": ",cosine_reconst))
  par(mfrow=c(length(alpha)+2,1))
  par(mar=c(1,2,4,1))
  barplot(hdp_sigs[,s], col=mut.cols, main=paste0("HDP ",s),names.arg="")
  barplot(reconsbs, col=mut.cols, main=paste0("Reconstituted ",s," cosine similarity to original: ", round(cosine_reconst, digits=2)))
  for (g in gdsigs) {
    barplot(ref[,g], col=mut.cols, main=paste0("PCAWG ", g, " accounts for ", round(alpha[g], digits=2)))
  }
}
dev.off()
# ref_sigs_R2 = ref_sigs_R2[ref_sigs_R2!='APOBEC']
# ref_sigs_R2 = unique(c(ref_sigs_R2,'SBS2','SBS13'))
# sigs_deconv_R2$Signature11 = c('SBS2','SBS13')
saveRDS(sigs_deconv_R2, file.path(HDP_DIR, "hdp2refsigs.Rdata"))
final_sigs=ref[,ref_sigs_R2]
write.table(final_sigs, file.path(HDP_DIR, "final_sigs.txt"), quote=F)
write.table(signature_fractionR2, file.path(HDP_DIR, "signature_fractionR2.txt"), quote=F)


### --------------- Refit reference signatures ----------------
library(sigfit)
library(RColorBrewer)
library(ape)
library(ggtree)
data("cosmic_signatures_v2")

final_sigs=t(read.table(file.path(HDP_DIR, "final_sigs.txt"), check.names=FALSE))

# --- SBSC -> SBS112 compatibility -------------------------------------------
# The refit LOGIC below is untouched; this only makes the names it works with
# agree with the renamed plot/summary labels. Both of its name sources are
# written by earlier runs and still say SBSC: final_sigs.txt (its columns) and
# hdp2refsigs.Rdata (its values). Without normalising both, line ~287's
# `sig_order[ref_sigs_present]` returns NA for SBSC and sort() drops it, so the
# signature would vanish from the refit with no error.
rownames(final_sigs)[rownames(final_sigs) == "SBSC"] <- "SBS112"
hdp_counts=read.table(SBS_INPUT, check.names=F)
hdp_counts=hdp_counts[apply(hdp_counts,1,sum)>50,]

hdp_exposures <- read.csv(file.path(HDP_DIR, "combined_HDP_exposure.csv"), row.names=1, check.names=F)
hdp_exposures <- hdp_exposures[,-1]   # Remove Component 0
hdp_exposures <- hdp_exposures[rownames(hdp_counts), , drop=FALSE]
colnames(hdp_exposures) <- paste0('Signature', 1:ncol(hdp_exposures))

hdp2ref=readRDS(file.path(HDP_DIR, "hdp2refsigs.Rdata"))
hdp2ref <- lapply(hdp2ref, function(x) { x[x == "SBSC"] <- "SBS112"; x })

# sig_order=1:10
# names(sig_order)=c("SBS1","SBS5","SBS18","SBS2","SBS13","SBS88",'SBS89','SBS112',"SBS17b",'SBS93')

sig_order=1:8
names(sig_order)=c("SBS1","SBS5","SBS18","SBS88",'SBS89','SBS112',"SBS17b",'SBS93')


getPalette = colorRampPalette(brewer.pal(9, "Set3"))
all_cols=getPalette(9)
all_cols=all_cols[c(1:3,6:8)]
all_cols=c(all_cols,"magenta","firebrick")
final_sigs = final_sigs[names(sort(sig_order[rownames(final_sigs)])),]
# names(all_cols)=c("SBS1","SBS5","SBS18","SBS2","SBS13","SBS88",'SBS89','SBS112',"SBS17b",'SBS93')
names(all_cols)=c("SBS1","SBS5","SBS18","SBS88",'SBS89','SBS112',"SBS17b",'SBS93')

patients=unique(sapply(strsplit(rownames(hdp_counts),'_'), function(x) x[1]))
colnames(hdp_counts)=colnames(cosmic_signatures_v2)

for (patient in patients[1:length(patients)]){
  hdp_exposures_tmp=hdp_exposures[rownames(hdp_counts)[rowSums(hdp_counts)>200], , drop=FALSE]
  hdp_sigs_final=colnames(hdp_exposures_tmp)[colSums(hdp_exposures_tmp[grepl(paste0(patient,"_"),rownames(hdp_exposures_tmp)),]>0.05)>0]
  for (sample in grep(paste0(patient,"_"),rownames(hdp_counts),value=T)){
    if (file.exists(file.path(HDP_DIR, "refit", paste0("fit_branch_",sample,".txt")))) next
    if (sample %in% rownames(hdp_exposures)){
    hdp_sigs_present=colnames(hdp_exposures)[hdp_exposures[sample,]>0.05]
    hdp_sigs_present=hdp_sigs_present[hdp_sigs_present %in% hdp_sigs_final]

    ref_sigs_present=unique(c(unlist(hdp2ref[hdp_sigs_present])))
    ref_sigs_present = ref_sigs_present[ref_sigs_present %in% rownames(final_sigs)]

    # consider SBS1/5/18 for every crypt
    # ref_sigs_present = unique(c(ref_sigs_present,'SBS1','SBS5','SBS18'))
    }
    else {
      ref_sigs_present = c('SBS1','SBS5','SBS18')
    }

    counts_patient=hdp_counts[sample,]
    colnames(counts_patient)=colnames(cosmic_signatures_v2)
    ref_sigs_present=names(sort(sig_order[ref_sigs_present]))

    fit=fit_signatures(counts=counts_patient,signatures = final_sigs[ref_sigs_present,],
                       iter = 20000,warmup = 10000,model="poisson",chains = 5)
    if (length(ref_sigs_present) > 1) {
      pars <- retrieve_pars(fit,par = "exposures",hpd_prob = 0.95)
      exposure_matrix<-t(pars$mean[,ref_sigs_present])
    } else {
      exposure_matrix <- data.frame(matrix(ncol=1, nrow=1, 1))
      colnames(exposure_matrix) = sample
      rownames(exposure_matrix) = ref_sigs_present
    }

    ref_sigs_present_final=rownames(exposure_matrix)[exposure_matrix[,]>0.05]

    # if (sum(exposure_matrix[rownames(exposure_matrix) %in% c('SBS2','SBS13'),])>0.05){
    #   ref_sigs_present_final=unique(c(ref_sigs_present_final,"SBS2","SBS13"))
    # }

    ref_sigs_present_final=names(sort(sig_order[ref_sigs_present_final]))
    if (length(ref_sigs_present_final) != length(ref_sigs_present)){
      fit=fit_signatures(counts=counts_patient,signatures = final_sigs[ref_sigs_present_final,],
                         iter = 20000,warmup = 10000,model="poisson",chains = 5)
      if (length(ref_sigs_present_final) > 1) {
        pars <- retrieve_pars(fit,par = "exposures",hpd_prob = 0.95)
        exposure_matrix<-t(pars$mean[,ref_sigs_present_final])
      } else {
        exposure_matrix <- data.frame(matrix(ncol=1, nrow=1, 1))
        colnames(exposure_matrix) = sample
        rownames(exposure_matrix) = ref_sigs_present_final
      }
    }
    write.table(exposure_matrix, file.path(HDP_DIR, "refit", paste0("fit_branch_",sample,".txt")), quote=F, sep='\t')

      try(plot_all(mcmc_samples = fit,
             out_path = file.path(HDP_DIR, "plots"),
             prefix = paste0("Fitting_",sample)))
  }
}

# -----------------Plot signatures onto phylogenetic trees-----------------
fap_patients  <- read.table(FAP_PATIENTS,  header=FALSE)$V1
norm_patients <- read.table(NORM_PATIENTS, header=FALSE)$V1

get_group <- function(sample_names, patient_list) {
  sample_names[sapply(strsplit(sample_names, '_'), function(x) x[1]) %in% patient_list]
}

# fap
hdp_counts=read.table(SBS_INPUT, check.names=F)
hdp_counts=hdp_counts[apply(hdp_counts,1,sum)>50,]
hdp_counts = hdp_counts[get_group(rownames(hdp_counts), fap_patients), , drop=FALSE]
patients    = fap_patients[fap_patients %in% unique(sapply(strsplit(rownames(hdp_counts),'_'), function(x) x[1]))]
colnames(hdp_counts)=colnames(cosmic_signatures_v2)
database = read.csv(DATABASE, stringsAsFactors=FALSE)

for (patient in patients){
  tree=try(read.tree(file.path(FAP_SNP_DIR, patient, paste0(patient,"_snp_tree_with_branch_length.tree"))))
  tree[["tip.label"]]=sapply(tree[["tip.label"]],function(x){
    label = database$group_id[database$sample==x]
    if (length(label)==0 || is.na(label) || label %in% c('normal','')) {
      label=""
    }
    label
  })
  tree_df=fortify(tree)

  samples = grep(paste0(patient,"_"),rownames(hdp_counts),value=T)
  exposure_matrix_sample=data.frame(matrix(ncol = length(samples), nrow = length(sig_order)))
  rownames(exposure_matrix_sample) <- names(sig_order)
  colnames(exposure_matrix_sample) <- samples

  for (sample in samples){
    exposures = read.table(file.path(HDP_DIR, "refit", paste0("fit_branch_",sample,".txt")), check.names=FALSE)
    # refit files written before the SBSC -> SBS112 rename still carry the old
    # row name; without this the `sig %in% rownames` test below falls through
    # to its `else 0` branch and SBS112 is silently zeroed everywhere
    rownames(exposures)[rownames(exposures) == "SBSC"] <- "SBS112"

    for (sig in names(sig_order)){
      if (sig %in% rownames(exposures)){
        exposure_matrix_sample[sig,sample] = exposures[sig,sample]
      }
      else{
        exposure_matrix_sample[sig,sample]=0
      }
    }
  }

  exposure_matrix_sample = exposure_matrix_sample[rowSums(exposure_matrix_sample)>0,]
  write.table(exposure_matrix_sample, file.path(HDP_DIR, "refit", paste0("fit_branch_",patient,".txt")), sep='\t', quote=F)
  cols=all_cols[rownames(exposure_matrix_sample)]
  ref_sigs_present_final=rownames(exposure_matrix_sample)
  branches=sapply(strsplit(samples,split="_"), function(x) x[[2]])

  n_tips <- length(tree$tip.label)
  pdf_h  <- max(7, n_tips * 0.25)
  pdf_w  <- max(9, pdf_h * 0.7)
  pdf(file.path(HDP_DIR, "tree_plots", "fap", paste0(patient,"_tree_with_sigs.pdf")), width=pdf_w, height=pdf_h)
  plot(tree,label.offset=0.01*max(tree_df$x),show.tip.label=T)
  for (k in 1:length(samples)){
    n=as.numeric(branches[k])
    x_end=tree_df$x[n]
    x_start=tree_df$x[tree_df$parent[n]]
    x_intv=x_end-x_start
    y=node.height(tree)[n]
    tipnum=sum(tree_df$isTip)
    for (s in ref_sigs_present_final){
      x_end=x_start+exposure_matrix_sample[s,samples[k]]*x_intv
      rect(ybottom=y-min(0.02*tipnum,0.4),ytop=y+min(0.02*tipnum,0.4),xleft=x_start,xright=x_end,col=cols[s],lwd=0.5)
      x_start=x_end
    }
  }
  axisPhylo(side = 1,backward=F)
  title(xlab=patient)
  legend("topright",title="Signatures", legend=ref_sigs_present_final, fill=cols, bty="n",cex=0.8, ncol=1, xjust=0.5)
  dev.off()
}


# control
hdp_counts=read.table(SBS_INPUT, check.names=F)
hdp_counts=hdp_counts[apply(hdp_counts,1,sum)>50,]
hdp_counts = hdp_counts[get_group(rownames(hdp_counts), norm_patients), , drop=FALSE]
patients    = norm_patients[norm_patients %in% unique(sapply(strsplit(rownames(hdp_counts),'_'), function(x) x[1]))]
colnames(hdp_counts)=colnames(cosmic_signatures_v2)

for (patient in patients){
  tree=try(read.tree(file.path(NORM_SNP_DIR, patient, paste0(patient,"_snp_tree_with_branch_length.tree"))))
  if (inherits(tree, "try-error")) { message("Tree not found: ", patient); next }
  tree_df=fortify(tree)

  samples = grep(paste0(patient,"_"),rownames(hdp_counts),value=T)
  exposure_matrix_sample=data.frame(matrix(ncol = length(samples), nrow = length(sig_order)))
  rownames(exposure_matrix_sample) <- names(sig_order)
  colnames(exposure_matrix_sample) <- samples

  for (sample in samples){
    exposures = read.table(file.path(HDP_DIR, "refit", paste0("fit_branch_",sample,".txt")), check.names=FALSE)
    # refit files written before the SBSC -> SBS112 rename still carry the old
    # row name; without this the `sig %in% rownames` test below falls through
    # to its `else 0` branch and SBS112 is silently zeroed everywhere
    rownames(exposures)[rownames(exposures) == "SBSC"] <- "SBS112"

    for (sig in names(sig_order)){
      if (sig %in% rownames(exposures)){
        exposure_matrix_sample[sig,sample] = exposures[sig,sample]
      }
      else{
        exposure_matrix_sample[sig,sample]=0
      }
    }
  }

  exposure_matrix_sample = exposure_matrix_sample[rowSums(exposure_matrix_sample)>0,]
  write.table(exposure_matrix_sample, file.path(HDP_DIR, "refit", paste0("fit_branch_",patient,"_control.txt")), sep='\t', quote=F)
  cols=all_cols[rownames(exposure_matrix_sample)]
  ref_sigs_present_final=rownames(exposure_matrix_sample)
  branches=sapply(strsplit(samples,split="_"), function(x) x[[2]])

  n_tips <- length(tree$tip.label)
  pdf_h  <- max(7, n_tips * 0.25)
  pdf_w  <- max(9, pdf_h * 0.7)
  pdf(file.path(HDP_DIR, "tree_plots", "control", paste0("control_",patient,"_tree_with_ref.pdf")), width=pdf_w, height=pdf_h)
  plot(tree,label.offset=0.01*max(tree_df$x),show.tip.label=T)
  for (k in 1:length(samples)){
    n=as.numeric(branches[k])
    x_end=tree_df$x[n]
    x_start=tree_df$x[tree_df$parent[n]]
    x_intv=x_end-x_start
    y=node.height(tree)[n]
    tipnum=sum(tree_df$isTip)
    for (s in ref_sigs_present_final){
      x_end=x_start+exposure_matrix_sample[s,samples[k]]*x_intv
      rect(ybottom=y-min(0.02*tipnum,0.4),ytop=y+min(0.02*tipnum,0.4),xleft=x_start,xright=x_end,col=cols[s],lwd=0.5)
      x_start=x_end
    }
  }
  axisPhylo(side = 1,backward=F)
  title(xlab=patient)
  legend("topright",title="Signatures", legend=ref_sigs_present_final, fill=cols, bty="n",cex=0.8, ncol=1, xjust=0.5)
  dev.off()
}

# -----------build exposure matrix-----------
sig_order=1:8
names(sig_order)=c("SBS1","SBS5","SBS18","SBS88",'SBS89','SBS112',"SBS17b",'SBS93')

# fap
hdp_counts=read.table(SBS_INPUT, check.names=F)
hdp_cohort1=hdp_counts[apply(hdp_counts,1,sum)>50,]
hdp_cohort1 = hdp_cohort1[get_group(rownames(hdp_cohort1), fap_patients), , drop=FALSE]
patients = fap_patients[fap_patients %in% unique(sapply(strsplit(rownames(hdp_cohort1),'_'), function(x) x[1]))]


exposure_matrix_all = data.frame(matrix(nrow = 0, ncol = length(sig_order)))
for (patient in patients){
  samples = grep(paste0(patient,"_"),rownames(hdp_cohort1),value=T)
  exposure_matrix_sample=data.frame(matrix(nrow = length(samples), ncol = length(sig_order)))
  colnames(exposure_matrix_sample) <- names(sig_order)
  rownames(exposure_matrix_sample) <- samples

  for (sample in samples){
    exposures = read.table(file.path(HDP_DIR, "refit", paste0("fit_branch_",sample,".txt")), check.names=FALSE)
    # refit files written before the SBSC -> SBS112 rename still carry the old
    # row name; without this the `sig %in% rownames` test below falls through
    # to its `else 0` branch and SBS112 is silently zeroed everywhere
    rownames(exposures)[rownames(exposures) == "SBSC"] <- "SBS112"

    for (sig in names(sig_order)){
      if (sig %in% rownames(exposures)){
        exposure_matrix_sample[sample,sig] = exposures[sig,sample]
      }
      else{
        exposure_matrix_sample[sample,sig]=0
      }
    }
  }
  exposure_matrix_all = rbind(exposure_matrix_all,exposure_matrix_sample)
}


# control
hdp_counts=read.table(SBS_INPUT, check.names=F)
hdp_cohort2=hdp_counts[apply(hdp_counts,1,sum)>50,]
hdp_cohort2 = hdp_cohort2[get_group(rownames(hdp_cohort2), norm_patients), , drop=FALSE]
patients = norm_patients[norm_patients %in% unique(sapply(strsplit(rownames(hdp_cohort2),'_'), function(x) x[1]))]

for (patient in patients){
  samples = grep(paste0(patient,"_"),rownames(hdp_cohort2),value=T)
  exposure_matrix_sample=data.frame(matrix(nrow = length(samples), ncol = length(sig_order)))
  colnames(exposure_matrix_sample) <- names(sig_order)
  rownames(exposure_matrix_sample) <- samples

  for (sample in samples){
    exposures = read.table(file.path(HDP_DIR, "refit", paste0("fit_branch_",sample,".txt")), check.names=FALSE)
    # refit files written before the SBSC -> SBS112 rename still carry the old
    # row name; without this the `sig %in% rownames` test below falls through
    # to its `else 0` branch and SBS112 is silently zeroed everywhere
    rownames(exposures)[rownames(exposures) == "SBSC"] <- "SBS112"

    for (sig in names(sig_order)){
      if (sig %in% rownames(exposures)){
        exposure_matrix_sample[sample,sig] = exposures[sig,sample]
      }
      else{
        exposure_matrix_sample[sample,sig]=0
      }
    }
  }
  exposure_matrix_all = rbind(exposure_matrix_all,exposure_matrix_sample)
}


write.table(exposure_matrix_all, file.path(HDP_DIR, "refit", "sig_exposure.txt"), sep='\t', quote=F)
