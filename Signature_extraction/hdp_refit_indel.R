library(RColorBrewer)
library(ggtree)
library(lsa)
library(lattice)
library(sigfit)
library(ape)

# --- Path configuration ---
HDP_DIR         <- "/lustre/scratch125/casm/teams/team267/users/yw2/colon/HDP/indel"
INDEL_INPUT     <- "/lustre/scratch125/casm/teams/team267/users/yw2/colon/LCMB_analysis/indel_fap_control_perbranch.hdp.txt"
FAP_INDEL_DIR   <- "/lustre/scratch125/casm/teams/team267/users/yw2/colon/LCMB_analysis/FAP/indel"
DATABASE        <- "/lustre/scratch125/casm/teams/team267/users/yw2/colon/LCMB_analysis/FAP/database_latest.txt"

dir.create(file.path(HDP_DIR, "refit"),                 showWarnings=FALSE, recursive=TRUE)
dir.create(file.path(HDP_DIR, "plots"),                 showWarnings=FALSE, recursive=TRUE)
dir.create(file.path(HDP_DIR, "tree_plots", "fap"),     showWarnings=FALSE, recursive=TRUE)
dir.create(file.path(HDP_DIR, "tree_plots", "control"), showWarnings=FALSE, recursive=TRUE)

# Load reference signatures
data("cosmic_signatures_indel_v3.2", package="sigfit")
ref <- t(cosmic_signatures_indel_v3.2)

mut.cols <- c(rep("#FDBE6F", 6), rep("#FF8001", 6), rep("#B0DD8B", 6), rep("#36A12E", 6),
              rep("#FDCAB5", 6), rep("#FC8A6A", 6), rep("#F14432", 6), rep("#BC141A", 6),
              rep("#D0E1F2", 6), rep("#94C4DF", 6), rep("#4A98C9", 6), rep("#1764AB", 6),
              rep("#E2E2EF", 1), rep("#B6B6D8", 2), rep("#8683BD", 3), rep("#61409B", 5))

# Load HDP signatures (TSV: rows = ID83 types, cols = X0-X5; header has n-1 names → auto row.names)
hdp_sigs <- read.table(file.path(HDP_DIR, "combined_HDP_sigs.txt"), check.names=F)
colnames(hdp_sigs) <- paste0('Signature', 0:(ncol(hdp_sigs)-1))
hdp_sigs <- hdp_sigs[,-1]   # Remove Signature0 (noise component)

# Load HDP exposures
hdp_exposures <- read.csv(file.path(HDP_DIR, "combined_HDP_exposure.csv"), row.names=1, check.names=F)
hdp_exposures <- hdp_exposures[,-1]   # Remove Component 0
colnames(hdp_exposures) <- paste0('Signature', 1:ncol(hdp_exposures))

input_for_hdp <- read.table(INDEL_INPUT, check.names=F)
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

write.table(cosine_matrix, file.path(HDP_DIR, "Cosine_similarities_pcawag.txt"), sep="\t", quote=F)

# Plot output
color.palette = colorRampPalette(c("white", "orange", "purple"))
pdf(file.path(HDP_DIR, "cosine_similarity.pdf"), width=9, height=4)
levelplot(t(cosine_matrix[dim(cosine_matrix)[1]:1,]),col.regions=color.palette, aspect="fill", scales=list(x=list(rot=90)))
dev.off()

# ----------Assess cosine similarities for all suspected sigs -------------
gdsigs=c("ID1","ID2","ID5","ID18","ID14")

ref=t(cosmic_signatures_indel_v3.2)
ref=cbind(ref[,gdsigs], hdp_sigs$Signature3, hdp_sigs$Signature5)
colnames(ref)[(ncol(ref)-1):ncol(ref)] = c('IDA','IDB')
gdsigs=unique(c(gdsigs,'IDA','IDB'))
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

# Plot output
color.palette = colorRampPalette(c("white", "orange", "purple"))
pdf(file.path(HDP_DIR, "cosine_similarity_ref.pdf"), width=9, height=4)
levelplot(t(cosine_matrix[dim(cosine_matrix)[1]:1,]),col.regions=color.palette, aspect="fill", scales=list(x=list(rot=90)))
dev.off()


# --------------First iteration; decomposed hdp sigs into all suspected signatures---------
cosine_matrix=data.frame(matrix(nrow=ncol(hdp_sigs), ncol=ncol(ref)))
rownames(cosine_matrix)=colnames(hdp_sigs)
colnames(cosine_matrix)=colnames(ref)

for (n in 1:nrow(cosine_matrix)) {
  for (m in 1:ncol(cosine_matrix)) {
    cosine_matrix[n,m] <- cosine(x=hdp_sigs[,rownames(cosine_matrix)[n]],
                                 y=ref[,colnames(cosine_matrix)[m]])
  }
}

sigs_to_decompose=rowSums(cosine_matrix>0.9)
for(n in 1:nrow(cosine_matrix)){
  print(paste0(rownames(cosine_matrix)[n],": ",paste(colnames(cosine_matrix)[cosine_matrix[n,]>0.9],collapse=",")))
}

signatures=t(ref[,gdsigs])

signature_fractionR1 = matrix(NA,nrow=nrow(signatures),ncol=length(colnames(hdp_sigs)))
rownames(signature_fractionR1) = rownames(signatures)
colnames(signature_fractionR1) = colnames(hdp_sigs)
maxiter <- 1000

for (j in 1:length(colnames(hdp_sigs))) {
  freqs = hdp_sigs[,j]
  freqs[is.na(freqs)] = 0
  alpha = runif(nrow(signatures)); alpha=alpha/sum(alpha)
  for (iter in 1:maxiter) {
    contr = t(array(alpha,dim=c(nrow(signatures),83))) * t(signatures)
    probs = contr/array(rowSums(contr),dim=dim(contr))
    probs = probs * freqs
    probs[is.na(probs)] = 0
    old_alpha = alpha
    alpha = colSums(probs)/sum(probs)
    if (sum(abs(alpha-old_alpha))<1e-5) {
      break
    }
  }
  signature_fractionR1[,j] = alpha
}

# Plot initial deconvolution
pdf(file.path(HDP_DIR, "sigs_decomposition_r1.pdf"), width=9, height=4)
color.palette = colorRampPalette(c("white", "orange", "purple"))
levelplot((signature_fractionR1[nrow(signature_fractionR1):1,]),col.regions=color.palette, aspect="fill", scales=list(x=list(rot=90)))
dev.off()

# Select which signatures to decompose into reference signatures
sigs_to_deconv=names(sigs_to_decompose)[sigs_to_decompose!=1]
sigs_to_deconv

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
for(s in colnames(hdp_sigs)){
  gdsigs <- sigs_deconv_R2[[s]]
  signatures <- t(ref[,gdsigs])

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

  alpha = runif(nrow(signatures)); alpha=alpha/sum(alpha)
  for (iter in 1:maxiter) {
    contr = t(array(alpha,dim=c(nrow(signatures),83))) * t(signatures)
    probs = contr/array(rowSums(contr),dim=dim(contr))
    probs = probs * freqs
    probs[is.na(probs)] = 0
    old_alpha = alpha
    alpha = colSums(probs)/sum(probs)
    if (sum(abs(alpha-old_alpha))<1e-5) {
      break
    }
  }
  signature_fractionR2[gdsigs,s] = alpha
  reconsbs <- rep(0,83)
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

saveRDS(sigs_deconv_R2, file.path(HDP_DIR, "hdp2refsigs.Rdata"))
final_sigs=ref[,ref_sigs_R2]
write.table(final_sigs, file.path(HDP_DIR, "final_sigs.txt"), quote=F)
write.table(signature_fractionR2, file.path(HDP_DIR, "signature_fractionR2.txt"), quote=F)


### --------------- Refit reference signatures ----------------
final_sigs=data.frame(t(read.table(file.path(HDP_DIR, "final_sigs.txt"), check.names=FALSE)), check.names=FALSE)
hdp_counts=read.table(INDEL_INPUT, check.names=F)
hdp_counts=hdp_counts[apply(hdp_counts,1,sum)>30,]

hdp_exposures <- read.csv(file.path(HDP_DIR, "combined_HDP_exposure.csv"), row.names=1, check.names=F)
hdp_exposures <- hdp_exposures[,-1]   # Remove Component 0
hdp_exposures <- hdp_exposures[rownames(hdp_counts), , drop=FALSE]
colnames(hdp_exposures) <- paste0('Signature', 1:ncol(hdp_exposures))

hdp2ref=readRDS(file.path(HDP_DIR, "hdp2refsigs.Rdata"))

sig_order=1:7
names(sig_order)=c("ID1","ID2","ID5","ID18","ID14","IDA","IDB")

final_sigs = final_sigs[names(sort(sig_order[rownames(final_sigs)])),]

patients=unique(sapply(strsplit(rownames(hdp_counts),'_'), function(x) x[1]))
for (patient in patients[1:length(patients)]){
  hdp_sigs_final=colnames(hdp_exposures)[colSums(hdp_exposures[grepl(paste0(patient,"_"),rownames(hdp_exposures)),]>0.05)>0]

  for (sample in grep(paste0(patient,"_"),rownames(hdp_counts),value=T)){
    if (file.exists(file.path(HDP_DIR, "refit", paste0("fit_",sample,".txt")))) next
    counts_patient=hdp_counts[sample,]
    if (sample %in% rownames(hdp_exposures)){
      # Only offer sigfit reference signatures that trace back to a raw HDP
      # component with genuine exposure (>0.05) in this branch -- do not
      # unconditionally force-add ID1/ID2/ID5 regardless of HDP support
      # (the original behaviour here injected spurious ID1/ID2/ID5 weight
      # into branches that HDP found to be ~ID18-only, e.g. OO10_7,
      # PD37449_19/35/37). Fall back to the default only when HDP found
      # nothing at all for this branch (handled in the else branch below).
      hdp_sigs_present=colnames(hdp_exposures)[hdp_exposures[sample,]>0.05]
      ref_sigs_present=unique(unlist(hdp2ref[hdp_sigs_present]))
      ref_sigs_present=ref_sigs_present[ref_sigs_present %in% rownames(final_sigs)]
      if (length(ref_sigs_present)==0) ref_sigs_present = c('ID1','ID2','ID5')
      ref_sigs_present=names(sort(sig_order[ref_sigs_present]))
    } else {
      ref_sigs_present = c('ID1','ID2','ID5')
    }

    fit=fit_signatures(counts=counts_patient, signatures=final_sigs[ref_sigs_present,],
                       iter=5000, warmup=2000, model="poisson", chains=5)

    if (length(ref_sigs_present)>1){
      pars <- retrieve_pars(fit, par="exposures", hpd_prob=0.95)
      exposure_matrix <- t(pars$mean[,ref_sigs_present])
    } else {
      exposure_matrix <- data.frame(matrix(ncol=1, nrow=1, 1))
      colnames(exposure_matrix) = sample
      rownames(exposure_matrix) = ref_sigs_present
    }

    ref_sigs_present_final=rownames(exposure_matrix)[exposure_matrix[,]>0.05]
    ref_sigs_present_final=names(sort(sig_order[ref_sigs_present_final]))

    if (!identical(ref_sigs_present_final, ref_sigs_present)){
      fit=fit_signatures(counts=counts_patient, signatures=final_sigs[ref_sigs_present_final,],
                         iter=5000, warmup=2000, model="poisson", chains=5)
      if (length(ref_sigs_present_final)>1){
        pars <- retrieve_pars(fit, par="exposures", hpd_prob=0.95)
        exposure_matrix <- t(pars$mean[,ref_sigs_present_final])
      } else {
        exposure_matrix <- data.frame(matrix(ncol=1, nrow=1, 1))
        colnames(exposure_matrix) = sample
        rownames(exposure_matrix) = ref_sigs_present_final
      }
    }
    write.table(exposure_matrix, file.path(HDP_DIR, "refit", paste0("fit_",sample,".txt")), quote=F, sep='\t')

    try(plot_all(mcmc_samples=fit,
                 out_path=file.path(HDP_DIR, "plots"),
                 prefix=paste0("Fitting_",sample)))
  }
}


# # -----------------Plot signatures onto phylogenetic trees-----------------
# sig_order=1:7
# names(sig_order)=c("ID1","ID2","ID18","ID14","IDA","IDB","IDC")

# all_cols=c("#DDDCC7","#FAEBEA","#DC978D","#F5B46F","#A6ADBA","#6993A4","#533A4D")
# names(all_cols)=c("ID1","ID2","ID18","ID14","IDA","IDB","IDC")

# # fap
# hdp_counts=read.table(INDEL_INPUT, check.names=F)
# hdp_cohort1 = hdp_counts[-grep('control', rownames(hdp_counts)),]
# patients=unique(sapply(strsplit(rownames(hdp_cohort1),'_'), function(x) x[1]))
# database = read.table(DATABASE, header=T, sep='\t')

# for (patient in patients){
#   tree=try(read.tree(file.path(FAP_INDEL_DIR, patient,
#                                paste0(patient,"_indel_on_snv_tree_with_branch_length.tree"))))
#   tree[["tip.label"]]=sapply(tree[["tip.label"]], function(x){
#     database$polyp_id[database$sample==x]
#   })
#   tree_df=fortify(tree)

#   samples = grep(paste0(patient,"_"), rownames(hdp_counts), value=T)
#   exposure_matrix_sample=data.frame(matrix(ncol=length(samples), nrow=length(sig_order)))
#   rownames(exposure_matrix_sample) <- names(sig_order)
#   colnames(exposure_matrix_sample) <- samples

#   for (sample in samples){
#     exposures = try(read.table(file.path(HDP_DIR, "refit", paste0("fit_",sample,".txt")), check.names=FALSE))
#     for (sig in names(sig_order)){
#       if (sig %in% rownames(exposures)){
#         exposure_matrix_sample[sig,sample] = exposures[sig,sample]
#       } else {
#         exposure_matrix_sample[sig,sample]=0
#       }
#     }
#   }
#   exposure_matrix_sample[is.na(exposure_matrix_sample)] = 0
#   exposure_matrix_sample = exposure_matrix_sample[rowSums(exposure_matrix_sample)>0,]
#   exposure_matrix_sample = exposure_matrix_sample[,colSums(exposure_matrix_sample)>0]
#   write.table(exposure_matrix_sample,
#               file.path(HDP_DIR, "refit", paste0("fit_",patient,".txt")), sep='\t', quote=F)

#   cols=all_cols[rownames(exposure_matrix_sample)]
#   ref_sigs_present_final=rownames(exposure_matrix_sample)
#   samples = samples[samples %in% colnames(exposure_matrix_sample)]
#   branches=sapply(strsplit(samples, split="_"), function(x) x[[2]])

#   pdf(file.path(HDP_DIR, "tree_plots", "fap", paste0(patient,"_tree_with_ref.pdf")), width=10, height=20)
#   plot(tree, label.offset=0.01*max(tree_df$x), show.tip.label=T)
#   for (k in 1:length(samples)){
#     n=as.numeric(branches[k])
#     x_end=tree_df$x[n]
#     x_start=tree_df$x[tree_df$parent[n]]
#     x_intv=x_end-x_start
#     y=node.height(tree)[n]
#     tipnum=sum(tree_df$isTip)
#     for (s in ref_sigs_present_final){
#       x_end=x_start+exposure_matrix_sample[s,samples[k]]*x_intv
#       rect(ybottom=y-min(0.02*tipnum,0.4), ytop=y+min(0.02*tipnum,0.4),
#            xleft=x_start, xright=x_end, col=cols[s], lwd=0.5)
#       x_start=x_end
#     }
#   }
#   axisPhylo(side=1, backward=F)
#   title(xlab=patient)
#   legend("topright", title="Signatures", legend=ref_sigs_present_final,
#          fill=cols, bty="n", cex=0.8, ncol=1, xjust=0.5)
#   dev.off()
# }


# -----------build exposure matrix-----------
sig_order=1:7
names(sig_order)=c("ID1","ID2","ID5","ID18","ID14","IDA","IDB")

hdp_counts=read.table(INDEL_INPUT, check.names=F)
hdp_cohort1 = hdp_counts[-grep('control', rownames(hdp_counts)),]
patients=unique(sapply(strsplit(rownames(hdp_cohort1),'_'), function(x) x[1]))

exposure_matrix_all = data.frame(matrix(nrow=0, ncol=length(sig_order)))
colnames(exposure_matrix_all) = names(sig_order)

for (patient in patients){
  samples = grep(paste0(patient,"_"), rownames(hdp_cohort1), value=T)
  exposure_matrix_sample=data.frame(matrix(nrow=length(samples), ncol=length(sig_order)))
  colnames(exposure_matrix_sample) <- names(sig_order)
  rownames(exposure_matrix_sample) <- samples

  for (sample in samples){
    exposures = try(read.table(file.path(HDP_DIR, "refit", paste0("fit_",sample,".txt")), check.names=FALSE))
    for (sig in names(sig_order)){
      if (sig %in% rownames(exposures)){
        exposure_matrix_sample[sample,sig] = exposures[sig,sample]
      } else {
        exposure_matrix_sample[sample,sig]=0
      }
    }
  }
  exposure_matrix_all = rbind(exposure_matrix_all, exposure_matrix_sample)
}

# control
hdp_cohort2 = hdp_counts[grep('control', rownames(hdp_counts)),]
patients=unique(sapply(strsplit(rownames(hdp_cohort2),'_'), function(x) x[1]))

for (patient in patients){
  samples = grep(paste0(patient,"_"), rownames(hdp_cohort2), value=T)
  exposure_matrix_sample=data.frame(matrix(nrow=length(samples), ncol=length(sig_order)))
  colnames(exposure_matrix_sample) <- names(sig_order)
  rownames(exposure_matrix_sample) <- samples

  for (sample in samples){
    exposures = try(read.table(file.path(HDP_DIR, "refit", paste0("fit_",sample,".txt")), check.names=FALSE))
    for (sig in names(sig_order)){
      if (sig %in% rownames(exposures)){
        exposure_matrix_sample[sample,sig] = exposures[sig,sample]
      } else {
        exposure_matrix_sample[sample,sig]=0
      }
    }
  }
  exposure_matrix_all = rbind(exposure_matrix_all, exposure_matrix_sample)
}

exposure_matrix_all[is.na(exposure_matrix_all)] = 0
exposure_matrix_all = exposure_matrix_all[apply(exposure_matrix_all,1,sum)>0,]
write.table(exposure_matrix_all, file.path(HDP_DIR, "refit", "sig_exposure.txt"), sep='\t', quote=F)
