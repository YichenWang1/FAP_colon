library(RColorBrewer)
library(lsa)
library(lattice)
library(sigfit)


# Load reference signatures
data("cosmic_signatures_v3.2")
features <- paste0(substr(colnames(cosmic_signatures_v3.2),1,1),'[',substr(colnames(cosmic_signatures_v3.2),2,2),'>',substr(colnames(cosmic_signatures_v3.2),6,6),']',substr(colnames(cosmic_signatures_v3.2),3,3))
ref<-t(cosmic_signatures_v3.2)
rownames(ref)<-features

mut.cols = rep(c("dodgerblue","black","red","grey70","olivedrab3","plum2"),each=16)

# Load HDP signatures
hdp_exposures=read.table("~/Documents/FARM/colon_hg38/snp/hdp/HDP_exposures.txt",check.names = F)
hdp_sigs=read.table("~/Documents/FARM/colon_hg38/snp/hdp/HDP_sigs.txt",check.names = F)
rownames(hdp_sigs)<-features

input_for_hdp<-read.table("~/Documents/FARM/colon_hg38/snp/hdp/sbs_fap_control_perbranch.txt",check.names = F)
input_for_hdp=input_for_hdp[apply(input_for_hdp,1,sum)>50,]

rownames(hdp_exposures)=rownames(input_for_hdp)
colnames(hdp_sigs) = paste0('Signature',c(0:(ncol(hdp_sigs)-1)))
rownames(hdp_sigs) = rownames(ref)
hdp_sigs=hdp_sigs[,-1]


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

write.table(cosine_matrix, "~/Documents/FARM/colon_hg38/snp/hdp/Cosine_similarities_cosmic.txt",sep="\t",quote=F)

# Plot output
color.palette = colorRampPalette(c("white", "orange", "purple"))
pdf("~/Documents/FARM/colon_hg38/snp/hdp/cosine_similarity.pdf",width=9, height=4)
levelplot(t(cosine_matrix[dim(cosine_matrix)[1]:1,]),col.regions=color.palette, aspect="fill", scales=list(x=list(rot=90)))
dev.off()

# ----------Assess cosine similarities for all suspected sigs -------------
# Use signatures in colon to complie a ist of possible reference signatures
gdsigs=c('SBS1', 'SBS2','SBS5', 'SBS13', 'SBS17a', 'SBS17b', 'SBS18', 'SBS28', 'SBS31', 'SBS35', 'SBS88', 'SBS89','SBS93')
# Add SBSC from Lee-Six 2019 into the signature list
Henry_sigs=read.table('../../../public/colon/Signature_extraction/SBS_novel_sig_trinucs.txt',header=T)

ref<-data.frame(t(cosmic_signatures_v3.2))
rownames(ref)<-features
ref=cbind(ref[,gdsigs],Henry_sigs$SBSC)
colnames(ref)[ncol(ref)] = 'SBSC'
gdsigs=c(gdsigs,"SBSC")
ref=ref[,gdsigs]
          
ref$APOBEC = (ref$SBS2 +ref$SBS13)/2

cosine_matrix=data.frame(matrix(nrow=ncol(hdp_sigs), ncol=ncol(ref)))
rownames(cosine_matrix)=colnames(hdp_sigs)
colnames(cosine_matrix)=colnames(ref)

for (n in 1:nrow(cosine_matrix)) {
  for (m in 1:ncol(cosine_matrix)) {
    cosine_matrix[n,m] <- cosine(x=hdp_sigs[,rownames(cosine_matrix)[n]],
                                 y=ref[,colnames(cosine_matrix)[m]])
  }
}

write.table(cosine_matrix, "~/Documents/FARM/colon_hg38/snp/hdp/Cosine_similarities_ref.txt",sep="\t",quote=F)


# Plot output
color.palette = colorRampPalette(c("white", "orange", "purple"))
pdf("~/Documents/FARM/colon_hg38/snp/hdp/cosine_similarity_ref.pdf",width=9, height=4)
levelplot(t(cosine_matrix[dim(cosine_matrix)[1]:1,]),col.regions=color.palette, aspect="fill", scales=list(x=list(rot=90)))
dev.off()


ref[,'SBSC']=hdp_sigs$Signature8

# --------------First iteration; decomposed hdp sigs into all suspected sigantures---------
sigs_to_decompose=rowSums(cosine_matrix>0.9)
for(n in 1:nrow(cosine_matrix)){
  print(paste0(rownames(cosine_matrix)[n],": ",paste(colnames(cosine_matrix)[cosine_matrix[n,]>0.9],collapse=",")))
}
# Select which signatures to decompose into reference signatures
sigs_to_deconv=names(sigs_to_decompose)[sigs_to_decompose!=1]
sigs_to_deconv

signatures=t(ref[,gdsigs])
signatures = signatures[rownames(signatures)!='SBSC',]  #Do not consider decomposition into SBSC
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
pdf("~/Documents/FARM/colon_hg38/snp/hdp/sigs_decomposition_r1.pdf",width=9, height=4)
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

pdf("~/Documents/FARM/colon_hg38/snp/hdp/HDPsig_decomposed.pdf",width=9, height=4)
# Repeat the deconvolution with the identified constitutive signatures
for(s in colnames(hdp_sigs)){
  gdsigs <- sigs_deconv_R2[[s]]
  signatures <- t(ref[,gdsigs])
  signatures = signatures[rownames(signatures)!='SBSC',]
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
ref_sigs_R2 = ref_sigs_R2[ref_sigs_R2!='APOBEC']
ref_sigs_R2 = unique(c(ref_sigs_R2,'SBS2','SBS13'))
sigs_deconv_R2$Signature11 = c('SBS2','SBS13')
saveRDS(sigs_deconv_R2,"~/Documents/FARM/colon_hg38/snp/hdp/hdp2refsigs.Rdata")
final_sigs=ref[,ref_sigs_R2]
write.table(final_sigs,"~/Documents/FARM/colon_hg38/snp/hdp/final_sigs.txt",quote = F)
write.table(signature_fractionR2,"~/Documents/FARM/colon_hg38/snp/hdp/signature_fractionR2.txt",quote = F)


### --------------- Refit reference signatures ----------------
library(sigfit)
library(RColorBrewer)
library(ape)
library(ggtree)
data("cosmic_signatures_v2")

final_sigs=t(read.table("~/Documents/FARM/colon_hg38/snp/hdp/final_sigs.txt",check.names=FALSE))
hdp_counts=read.table("~/Documents/FARM/colon_hg38/snp/hdp/sbs_fap_control_perbranch.txt",check.names = F)
hdp_counts=hdp_counts[apply(hdp_counts,1,sum)>50,]

hdp_exposures=read.table("~/Documents/FARM/colon_hg38/snp/hdp/HDP_exposures.txt",check.names = F)
hdp_exposures=hdp_exposures[,-1] # Remove Component 0
rownames(hdp_exposures)=rownames(hdp_counts)
colnames(hdp_exposures)= paste0('Signature',c(1:(ncol(hdp_exposures))))


hdp2ref=readRDS("~/Documents/FARM/colon_hg38/snp/hdp/hdp2refsigs.Rdata")

sig_order=1:10
names(sig_order)=c("SBS1","SBS5","SBS18","SBS2","SBS13","SBS88",'SBS89','SBSC',"SBS17b",'SBS93')

getPalette = colorRampPalette(brewer.pal(9, "Set3"))
all_cols=getPalette(9)
all_cols=all_cols[1:8]
all_cols=c(all_cols,"magenta","firebrick")
final_sigs = final_sigs[names(sort(sig_order[rownames(final_sigs)])),]
names(all_cols)=c("SBS1","SBS5","SBS18","SBS2","SBS13","SBS88",'SBS89','SBSC',"SBS17b",'SBS93')

# all_cols

patients=unique(sapply(strsplit(rownames(hdp_counts),'_'), function(x) x[1]))
colnames(hdp_counts)=colnames(cosmic_signatures_v2)

for (patient in patients[1:length(patients)]){
  hdp_exposures_tmp=hdp_exposures[rowSums(hdp_counts)>200,]
  hdp_sigs_final=colnames(hdp_exposures_tmp)[colSums(hdp_exposures_tmp[grepl(paste0(patient,"_"),rownames(hdp_exposures_tmp)),]>0.05)>0]
  for (sample in grep(paste0(patient,"_"),rownames(hdp_counts),value=T)){
    if (sample %in% rownames(hdp_exposures)){
    hdp_sigs_present=colnames(hdp_exposures)[hdp_exposures[sample,]>0.05]
    hdp_sigs_present=hdp_sigs_present[hdp_sigs_present %in% hdp_sigs_final]
    
    ref_sigs_present=unique(c(unlist(hdp2ref[hdp_sigs_present])))
    ref_sigs_present = ref_sigs_present[ref_sigs_present %in% rownames(final_sigs)]
    
    # consider SBS1/5/18 for every crypts
    ref_sigs_present = unique(c(ref_sigs_present,'SBS1','SBS5','SBS18'))
    }
    else {
      ref_sigs_present = c('SBS1','SBS5','SBS18')
    }
    
    counts_patient=hdp_counts[sample,]
    colnames(counts_patient)=colnames(cosmic_signatures_v2)
    ref_sigs_present=names(sort(sig_order[ref_sigs_present]))
    
    fit=fit_signatures(counts=counts_patient,signatures = final_sigs[ref_sigs_present,],
                       iter = 20000,warmup = 10000,model="poisson",chains = 5)
    pars <- retrieve_pars(fit,par = "exposures",hpd_prob = 0.95)
    exposure_matrix<-t(pars$mean[,ref_sigs_present])
    
    ref_sigs_present_final=rownames(exposure_matrix)[exposure_matrix[,]>0.05]
    
    if (sum(exposure_matrix[rownames(exposure_matrix) %in% c('SBS2','SBS13'),])>0.05){
      ref_sigs_present_final=unique(c(ref_sigs_present_final,"SBS2","SBS13"))
    }
    
    ref_sigs_present_final=names(sort(sig_order[ref_sigs_present_final]))
    if (length(ref_sigs_present_final) != length(ref_sigs_present)){
      fit=fit_signatures(counts=counts_patient,signatures = final_sigs[ref_sigs_present_final,],
                         iter = 20000,warmup = 10000,model="poisson",chains = 5)
      pars <- retrieve_pars(fit,par = "exposures",hpd_prob = 0.95)
      exposure_matrix<-t(pars$mean[,ref_sigs_present_final])
    }
    write.table(exposure_matrix,paste0("~/Documents/FARM/colon_hg38/snp/hdp/refit/fit_branch_",sample,".txt"),quote = F,sep='\t')
  
      plot_all(mcmc_samples = fit,
             out_path = "~/Documents/FARM/colon_hg38/snp/hdp/plots",
             prefix = paste0("Fitting_",sample))
  }
}

# -----------build exposure matrix-----------
sig_order=1:10
names(sig_order)=c("SBS1","SBS5","SBS18","SBS2","SBS13","SBS88",'SBS89','SBSC',"SBS17b",'SBS93')

# fap
hdp_counts=read.table("~/Documents/FARM/colon_hg38/snp/hdp/sbs_fap_control_perbranch.txt",check.names = F)
hdp_cohort1=hdp_counts[apply(hdp_counts,1,sum)>50,]
hdp_cohort1 = hdp_cohort1[-grep('control', rownames(hdp_cohort1)),]
patients=unique(sapply(strsplit(rownames(hdp_cohort1),'_'), function(x) x[1]))


exposure_matrix_all = data.frame(matrix(nrow = 0, ncol = length(sig_order)))
for (patient in patients){
  samples = grep(paste0(patient,"_"),rownames(hdp_cohort1),value=T)
  exposure_matrix_sample=data.frame(matrix(nrow = length(samples), ncol = length(sig_order)))
  colnames(exposure_matrix_sample) <- names(sig_order)
  rownames(exposure_matrix_sample) <- samples
  
  for (sample in samples){
    exposures = read.table(paste0("~/Documents/FARM/colon_hg38/snp/hdp/refit/fit_branch_",sample,".txt"),check.names=FALSE)
    
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
hdp_counts=read.table("~/Documents/FARM/colon_hg38/snp/hdp/sbs_fap_control_perbranch.txt",check.names = F)
hdp_cohort2=hdp_counts[apply(hdp_counts,1,sum)>50,]
hdp_cohort2 = hdp_cohort2[grep('control', rownames(hdp_cohort2)),]
patients=unique(sapply(strsplit(rownames(hdp_cohort2),'_'), function(x) x[1]))

for (patient in patients){
  samples = grep(paste0(patient,"_"),rownames(hdp_cohort2),value=T)
  exposure_matrix_sample=data.frame(matrix(nrow = length(samples), ncol = length(sig_order)))
  colnames(exposure_matrix_sample) <- names(sig_order)
  rownames(exposure_matrix_sample) <- samples
  
  for (sample in samples){
    exposures = read.table(paste0("~/Documents/FARM/colon_hg38/snp/hdp/refit/fit_branch_",sample,".txt"),check.names=FALSE)
    
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


write.table(exposure_matrix_all,"~/Documents/FARM/colon_hg38/snp/hdp/refit/sig_exposure.txt",sep = '\t',quote=F)

