library(ggtree)
library(ape)

# -----------------Plot signatures onto phylogenetic trees-----------------
# fap
hdp_counts=read.table("~/Documents/FARM/colon_hg38/snp/hdp/sbs_fap_control_perbranch.txt",check.names = F)
hdp_counts=hdp_counts[apply(hdp_counts,1,sum)>50,]
hdp_counts = hdp_counts[-grep('control', rownames(hdp_counts)),]
patients=unique(sapply(strsplit(rownames(hdp_counts),'_'), function(x) x[1]))
colnames(hdp_counts)=colnames(cosmic_signatures_v2)
database = read.csv('../data/fap_database_withexposure.csv',header=T)
for (patient in patients){
  tree=try(read.tree(paste0("~/Documents/FARM/colon_hg38/fap/snp/",patient,"/",patient,"_snp_tree_with_branch_length.tree")))
  tree[["tip.label"]]=sapply(tree[["tip.label"]],function(x){
    label = database$label_id[database$sample==x]
    if (label=='normal'){
      label=""
    }
    else {
      label=substr(label,8, nchar(label))
    }
    label
  })
  tree_df=fortify(tree)
  
  samples = grep(paste0(patient,"_"),rownames(hdp_counts),value=T)
  exposure_matrix_sample=data.frame(matrix(ncol = length(samples), nrow = length(sig_order)))
  rownames(exposure_matrix_sample) <- names(sig_order)
  colnames(exposure_matrix_sample) <- samples
  
  for (sample in samples){
    exposures = read.table(paste0("~/Documents/FARM/colon_hg38/snp/hdp/refit/fit_branch_",sample,".txt"),check.names=FALSE)
    
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
  write.table(exposure_matrix_sample,paste0("~/Documents/FARM/colon_hg38/snp/hdp/refit/fit_branch_",patient,".txt"),sep = '\t',quote=F,)
  cols=all_cols[rownames(exposure_matrix_sample)]
  ref_sigs_present_final=rownames(exposure_matrix_sample)
  branches=sapply(strsplit(samples,split="_"), function(x) x[[2]])
  
  pdf(paste0("~/Documents/FARM/colon_hg38/snp/tree_plots/fap/labelled/",patient,"_tree_with_sigs.pdf"),width = 7, height = 9) #h=7,9,18,w=7
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
hdp_counts=read.table("~/Documents/FARM/colon_hg38/snp/hdp/sbs_fap_control_perbranch.txt",check.names = F)
hdp_counts=hdp_counts[apply(hdp_counts,1,sum)>50,]
hdp_counts = hdp_counts[grep('control', rownames(hdp_counts)),]
patients=unique(sapply(strsplit(rownames(hdp_counts),'_'), function(x) x[1]))
colnames(hdp_counts)=colnames(cosmic_signatures_v2)

for (patient in patients[1:length(patients)]){
  # if (patient =='patient47') next
  
  tree=try(read.tree(paste0("~/Documents/FARM/colon_hg38/normal/snp/",patient,"/",patient,"_snp_tree_with_branch_length.tree")))
  tree_df=fortify(tree)
  samples = grep(paste0(patient,"_"),rownames(hdp_counts),value=T)
  exposure_matrix_sample=data.frame(matrix(ncol = length(samples), nrow = length(sig_order)))
  rownames(exposure_matrix_sample) <- names(sig_order)
  colnames(exposure_matrix_sample) <- samples
  
  for (sample in samples){
    exposures = read.table(paste0("~/Documents/FARM/colon_hg38/snp/hdp/refit/fit_branch_",sample,".txt"),check.names=FALSE)
    
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
  write.table(exposure_matrix_sample,paste0("~/Documents/FARM/colon_hg38/snp/hdp/refit/fit_branch_",patient,"_control.txt"),sep = '\t',quote=F,)
  cols=all_cols[rownames(exposure_matrix_sample)]
  ref_sigs_present_final=rownames(exposure_matrix_sample)
  branches=sapply(strsplit(samples,split="_"), function(x) x[[2]])
  
  pdf(paste0("~/Documents/FARM/colon_hg38/snp/tree_plots/control/control_",patient,"_tree_with_ref.pdf"),width = 7, height = 7)
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
