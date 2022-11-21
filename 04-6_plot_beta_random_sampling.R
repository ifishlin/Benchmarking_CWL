analysis_name <- "04-6_plot_beta_random_sampling"
config = file.path(getwd(), "00_project_setting.R")
source(config)

library(RnBeads)
library(ComplexUpset)
library(liftOver)
library(rtracklayer)
library(GenomicRanges)
library(dplyr)

fdr_th <- 0.1
for(prot in protocols_int){
  dmls <- read_(paste0("dml_i450_test_", prot), "04_DSS")
  dmrs <- read_(paste0("dmr_i450_list_", prot), "04_DSS")
  diffMethTable_site_cmp1 <- read.csv(file.path(data_dir, "diffMethTable_site_cmp1.csv"))
  
  #prepare GRanges from everything
  diffMethTable_site_cmp1$end <- diffMethTable_site_cmp1$Start+1
  diffMeth <- makeGRangesFromDataFrame(diffMethTable_site_cmp1, keep.extra.columns = T)
  
  for (dml in 1:length(dmls)){
    dmls[[dml]]$end <- dmls[[dml]]$pos+1
    colnames(dmls[[dml]])[2] <- "start"}
  
  dmls2 <- lapply(dmls, makeGRangesFromDataFrame, keep.extra.columns = T)
  
  #lift over the Illumina probes to hg38
  path = system.file(package="liftOver", "extdata", "hg19ToHg38.over.chain")
  ch = import.chain(path)
  seqlevelsStyle(diffMeth) = "UCSC"  # necessary
  cur38 = liftOver(diffMeth, ch)
  cur38 <-subset(cur38, lengths(cur38)==1)
  cur38 <- unlist(cur38)
  
  # extract all CpG sites
  #cpgs38 <- methrix::extract_CPGs(ref_genome = "BSgenome.Hsapiens.UCSC.hg38")
  #cpgs38 <- makeGRangesFromDataFrame(as.data.frame(cpgs38$cpgs))
  
  
  #prepare the data for the upset plot
  upst <- makeGRangesFromDataFrame(as.data.frame(cur38))
  down <- up <- upst
  up$Direction <- "Hypermethylated"
  down$Direction <- "Hypomethylated"
  up$IlluminaArray <- cur38$diffmeth.p.adj.fdr<fdr_th & cur38$mean.diff<0
  down$IlluminaArray <- cur38$diffmeth.p.adj.fdr<fdr_th & cur38$mean.diff>0
  
  
  for (dml in 1:length(dmls)){
    up$meth <- FALSE
    up$meth[queryHits(findOverlaps(up, dmls2[[dml]][dmls2[[dml]]$fdr<fdr_th & dmls2[[dml]]$diff<0,], type="equal"))] <- TRUE
    colnames(mcols(up))[colnames(mcols(up))=="meth"] <- gsub(paste0("_", replace_prot_name_(prot)), "", names(dmls)[dml])
    
    down$meth <- FALSE
    down$meth[queryHits(findOverlaps(down, dmls2[[dml]][dmls2[[dml]]$fdr<fdr_th & dmls2[[dml]]$diff>0,], type="equal"))] <- TRUE
    colnames(mcols(down))[colnames(mcols(down))=="meth"] <- gsub(paste0("_", replace_prot_name_(prot)), "", names(dmls)[dml])
  }
  
  ## add p value distribution for the first category for the workflows and compare it to a random set of sites.
  upst <- rbind(as.data.frame(up), as.data.frame(down))
  upst <- as.data.frame(upst)
  
  overlapping_sites <- upst[rowSums(upst[,-(1:6)])== ncol(upst[,-(1:6)]),] %>%
    distinct() %>%
    as.data.frame() %>%
    makeGRangesFromDataFrame(keep.extra.columns = T)
  
  df <- data.frame(pval=numeric(),
                   workflow=character())
  random <- data.frame(pval=numeric(),
                       workflow=character())
  
  for (dml in 1:length(dmls)){
    p.values_overlapping <- mcols(subsetByOverlaps(dmls2[[dml]], overlapping_sites))$pval
    p.values_random <- dmls2[[dml]]$pval[sample(x = which(dmls2[[dml]]$fdr<fdr_th), size = length(overlapping_sites), replace = F)]
    df <- rbind(df, data.frame(pval=p.values_overlapping, workflow= gsub(paste0("_", replace_prot_name_(prot)), "", names(dmls)[dml])))
    random <- rbind(random, data.frame(pval=p.values_random, workflow= gsub(paste0("_", replace_prot_name_(prot)), "", names(dmls)[dml])))
    #random <- rbind(random, data.frame(pval=p.values_random, workflow= names(dmls2)[[dml]]))
  }
  
  df$source <- "Shared"
  random$source <- "Random"
  df <- rbind(df, random)
  
  df$workflow = replace_wf_name_(df$workflow)
  
  df$pval2  <- -log10(df$pval)
  
  p <- ggplot(df, aes(x=pval2)) +
    geom_density(aes(color=source), size=1.3) + xlim(0,30) + facet_wrap(vars(workflow), nrow=2) + #theme_bw()
    #geom_density(aes(color=source), size=0.3) + xlim(0,30) + facet_wrap(vars(workflow), nrow=2) + #theme_bw()
    theme(panel.background = element_rect(fill = "white", colour = "black", linetype="solid"),
          panel.border=element_rect(fill=NA),
          strip.background = element_blank(),
          strip.text = element_text(size = 9),
          axis.text = element_text(size = 9),
          axis.title = element_text(size = 9), 
          #strip.text = element_text(size = 6),
          #axis.text = element_text(size = 5),
          #axis.title = element_text(size = 6),
          legend.position='right') + xlab("p-value")
  
  p
  
  save_(paste0("p_dist_incl_illumina_", prot), plot=p, use_pdf=TRUE,
        width=7,
        height=3, do_ggsave=FALSE) 

}