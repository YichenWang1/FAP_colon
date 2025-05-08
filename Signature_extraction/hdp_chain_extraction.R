#extract results from 10 HDP chains

library(hdp)
library(ggplot2)
library(sigfit)
input_for_hdp<-read.table("~/Documents/FARM/colon_hg38/snp/hdp/sbs_fap_control_perbranch.txt",check.names = F)
input_for_hdp=input_for_hdp[apply(input_for_hdp,1,sum)>50,]


chlist <- vector("list", 10) 
for (i in 1:10){
  chlist[[i]] <- readRDS(paste0("~/Documents/FARM/colon_hg38/snp/hdp/chlist_20241015_denovo_combined_",i,".rds"))

  }


mut_example_multi <- hdp_multi_chain(chlist)

lapply(chains(mut_example_multi), plot_lik, bty="L", start = 1000)
lapply(chains(mut_example_multi), plot_numcluster, bty="L")
lapply(chains(mut_example_multi), plot_data_assigned, bty="L")

mut_example_multi <- hdp_extract_components(mut_example_multi)
plot_comp_size(mut_example_multi, bty="L")

sub_vec = c("C>A","C>G","C>T","T>A","T>C","T>G")
ctx_vec = paste(rep(c("A","C","G","T"),each=4),rep(c("A","C","G","T"),times=4),sep="-")
full_vec = paste(rep(sub_vec,each=16),rep(ctx_vec,times=6),sep=",")
trinuc_context <- full_vec
group_factor <- as.factor(rep(c("C>A", "C>G", "C>T", "T>A", "T>C", "T>G"),
                              each=16))

mut_colours <- c(RColorBrewer::brewer.pal(10, 'Paired')[seq(1,10,2)], 'grey70')

pdf("~/Documents/FARM/colon_hg38/snp/hdp/HDP_sigs.pdf",width=8, height=4)
plot_comp_distn(mut_example_multi, cat_names=trinuc_context,
                grouping=group_factor, col=mut_colours,
                col_nonsig="grey80", show_group_labels=TRUE)# dev.off()
dev.off()

data("cosmic_signatures_v3.2")
features <- paste0(substr(colnames(cosmic_signatures_v3.2),1,1),'[',substr(colnames(cosmic_signatures_v3.2),2,2),'>',substr(colnames(cosmic_signatures_v3.2),6,6),']',substr(colnames(cosmic_signatures_v3.2),3,3))


hdp_exposures=mut_example_multi@comp_dp_distn[["mean"]][2:dim(mut_example_multi@comp_dp_distn[["mean"]])[1],]
rownames(hdp_exposures)=rownames(input_for_hdp)
write.table(hdp_exposures,"~/Documents/FARM/colon_hg38/snp/hdp/HDP_exposures.txt",quote = F)

hdp_sigs=t(mut_example_multi@comp_categ_distn[["mean"]][1:dim(mut_example_multi@comp_categ_distn[["mean"]])[1],])
rownames(hdp_sigs)<-features
write.table(hdp_sigs,"~/Documents/FARM/colon_hg38/snp/hdp/HDP_sigs.txt",quote = F)

