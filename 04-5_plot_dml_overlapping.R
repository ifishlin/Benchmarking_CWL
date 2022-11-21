analysis_name <- "04-5_plot_dml_overlapping"
config = file.path(getwd(), "0_project_setting.R")
source(config)

library(RnBeads)
library(ComplexUpset)
library(liftOver)
library(GenomicRanges)
library(rtracklayer)

fdr_th <- 0.1
for(prot in protocols_int) {
  print(prot)
  dmls <- read_(paste0("dml_i450_test_", prot), "04_DSS")
  dmrs <- read_(paste0("dmr_i450_list_", prot), "04_DSS")
  diffMethTable_site_cmp1 <- read.csv(file.path(data_dir, "diffMethTable_site_cmp1.csv"))
  
  #prepare GRanges from everything
  diffMethTable_site_cmp1$end <- diffMethTable_site_cmp1$Start+1
  diffMeth <- makeGRangesFromDataFrame(diffMethTable_site_cmp1, keep.extra.columns = T)
  
  
  for (dml in 1:length(dmls)){
    dmls[[dml]]$end <- dmls[[dml]]$pos+1
    colnames(dmls[[dml]])[2] <- "start"
  }
  
  dmls2 <- lapply(dmls, makeGRangesFromDataFrame, keep.extra.columns = T)  
  
  #lift over the Illumina probes to hg38
  path = system.file(package="liftOver", "extdata", "hg19ToHg38.over.chain")
  ch = import.chain(path)
  seqlevelsStyle(diffMeth) = "UCSC"  # necessary
  cur38 = liftOver(diffMeth, ch)
  cur38 <-subset(cur38, lengths(cur38)==1)
  cur38_tb = data.frame(unique(unlist(cur38))) #460499
  
  
  # extract all CpG sites
  #cpgs38 <- methrix::extract_CPGs(ref_genome = "BSgenome.Hsapiens.UCSC.hg38")
  #cpgs38 <- makeGRangesFromDataFrame(as.data.frame(cpgs38$cpgs))
  
  
  #prepare the data for the upset plot
  upst <- makeGRangesFromDataFrame(as.data.frame(cur38))
  down <- up <- upst
  up$Direction <- "Hypermethylated"
  down$Direction <- "Hypomethylated"
  a = data.frame(cur38)
  up$IlluminaArray <- a$diffmeth.p.adj.fdr<fdr_th & a$mean.diff<0 # 30380
  down$IlluminaArray <- a$diffmeth.p.adj.fdr<fdr_th & a$mean.diff>0 # 33846
  
  for (dml in 1:length(dmls)){
    up$meth <- FALSE
    up$meth[queryHits(findOverlaps(up, dmls2[[dml]][dmls2[[dml]]$fdr<fdr_th & dmls2[[dml]]$diff<0,], type="equal"))] <- TRUE
    colnames(mcols(up))[colnames(mcols(up))=="meth"] <- gsub(paste0("_", replace_prot_name_(prot)), "", names(dmls)[dml])
    
    down$meth <- FALSE
    down$meth[queryHits(findOverlaps(down, dmls2[[dml]][dmls2[[dml]]$fdr<fdr_th & dmls2[[dml]]$diff>0,], type="equal"))] <- TRUE
    colnames(mcols(down))[colnames(mcols(down))=="meth"] <- gsub(paste0("_", replace_prot_name_(prot)), "", names(dmls)[dml])
  }
  
  upst <- rbind(as.data.frame(up), as.data.frame(down))
  upst <- as.data.frame(upst)
  upst <- upst[rowSums(upst[,-(1:6)])>0,]
  plot = upset(
    upst,
    colnames(upst)[-(1:6)],
    # base_annotations=list(
    #   'Intersection size'=intersection_size(
    #     counts=FALSE,
    #     mapping=aes(fill=Direction))+
    #     scale_fill_lancet()
    # ),
    base_annotations=list(
      'Intersection size'=intersection_size(
        counts=TRUE,
        mapping=aes(fill=Direction),
        text=list(
          #vjust=-0.1,
          hjust=-0.1,
          angle=45,
          size=3,
          colour = "black"
        )
      )  + theme(legend.position='none') + color_palette_fill() + color_palette_color() #
    ),  
    # width_ratio=0.2, 
    # height_ratio=0.5, 
    # n_intersections=30, 
    # sort_sets="ascending",
    width_ratio=0.2,
    n_intersections=30,
    height_ratio=0.8,  
    matrix = intersection_matrix(
      geom=geom_point(
        size=1.2
      )
    ),
    set_sizes=(
      upset_set_size(
        geom=geom_bar(
          aes(fill=Direction),
          stat='count',
          width=0.6
        )
      )
      #+ theme(legend.position='none', axis.ticks.x = element_line()) 
      #+ scale_y_reverse(breaks=c(0, 1000, 2000, 3000, 4000),  label = c("0", "1k", "2k", "3k", "4k"))
      + color_palette_fill() + color_palette_color() #+ scale_y_break(c(34500, 65000)) + scale_y_break(c(3000, 33500))  
      #+ scale_y_reverse()
      + theme(legend.position='none', axis.text.x = element_text(size=8, angle = 60, vjust = 0.7), axis.ticks.x = element_line()) 
    ),
  )
  
  plot
  
  save_(paste0("DML_overlap_complexupset_IL_Seq_", prot), plot=plot, use_pdf=TRUE,
        width=7.5,
        height=4, do_ggsave=FALSE)   
}