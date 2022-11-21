library(tidyverse)
library(GenomicRanges)
library(ggplot2)
library(ggsci)

setwd("/home/y306n/OE0219YUYU/benchmarking/analysis/data_processing_manuscript")
proj_dir <- getwd()
figures_dir <- file.path(getwd(), "figures")
data_dir <- file.path(getwd(), "data")
if (!file.exists(figures_dir)){
  dir.create(figures_dir)
}
if (!file.exists(data_dir)){
  dir.create(data_dir)
}

figures_dir_ <- file.path(figures_dir, analysis_name)
data_dir_ <- file.path(data_dir, analysis_name)
if (!file.exists(figures_dir_)){
  dir.create(figures_dir_)
}
if (!file.exists(data_dir_)){
  dir.create(data_dir_)
}

dirname_2_intname <- function(protocol_name){
  protocol_name <- gsub("SWIFT", "Swift", protocol_name)
  protocol_name <- gsub("EMSEQ", "EM-seq", protocol_name)  
  protocol_name <- gsub("TWGBS", "T-WGBS", protocol_name)  
  return(protocol_name)
}


replace_prot_name_ <- function(protocol_name){
  protocol_name <- gsub("SWIFT", "Swift", protocol_name)
  protocol_name <- gsub("TWGBS", "T-WGBS", protocol_name)  
  protocol_name <- gsub("EMSEQ", "EM-seq", protocol_name)  
  protocol_name <- gsub("EMseq", "EM-seq", protocol_name)  
  return(protocol_name)
}

replace_wf_name_ <- function(workflow_name){
  # workflow_name <- gsub("Bismarkkersten", "Bismark", workflow_name)
  # workflow_name <- gsub("bwamethmbias", "bwa-meth", workflow_name)  
  # workflow_name <- gsub("Bismark_kersten_thesis", "Bismark", workflow_name)
  # workflow_name <- gsub("bwa-meth-mbias", "bwa-meth", workflow_name)
  workflow_name <- gsub("bwameth", "bwa-meth", workflow_name)
  return(workflow_name)
}

protocols_int <- c(
  "WGBS",
  "SWIFT",
  "TWGBS",
  "PBAT",
  "EMSEQ"
)

protocols <- c(
  "WGBS",
  "Swift",
  "T-WGBS",
  "PBAT",
  "EM-seq"
)

## methrix forbids "-" in the ID.
workflows_methrix <- c(
  "BAT",
  "Biscuit",
  "Bismark",
  "BSBolt",
  "bwameth",
  "FAME",
  "gemBS",
  "GSNAP",
  "methylCtools",
  "methylpy"
)

workflows_int = workflows_methrix

workflows <- c(
  "BAT",
  "Biscuit",
  "Bismark",
  "BSBolt",
  "bwa-meth",
  "FAME",
  "gemBS",
  "GSNAP",
  "methylCtools",
  "methylpy"
)

set_protocol_order <- function(tb){
  tb %>%
    dplyr::ungroup() %>%
    dplyr::mutate(protocol=ordered(protocol, levels=protocols))
}

samples <- c(
  "5N",
  "6N",
  "5T",
  "6T"
)

set_protocol_order <- function(tb){
  tb %>%
    dplyr::ungroup() %>%
    dplyr::mutate(protocol=ordered(protocol, levels=protocols_int))
}

set_sample_order <- function(tb){
  tb %>%
    dplyr::ungroup() %>%
    dplyr::mutate(sample=ordered(sample, levels=samples))
}

set_order <- function(tb){
  tb %>%
    set_protocol_order() %>%
    set_sample_order()
}

sample_color_sheme <- c(
  "#4592ff",
  "#fc7c38",
  "#1921ff",
  "#ff3700"
)
names(sample_color_sheme) <- samples


sample_linetype_sheme <- c(
  "dotted",
  "dotted",
  "dashed",
  "dashed"
)
names(sample_linetype_sheme) <- samples


sample_shape_sheme <- c(
  15,
  22,
  17,
  24
)
names(sample_shape_sheme) <- samples

sample_alpha_sheme <- c(
  0.5,
  1,
  0.5,
  1
)
names(sample_alpha_sheme) <- samples

protocol_color_sheme <- c(
  "#CD534CFF", #EB534C
  "#33cc00",
  "#0073C2FF", #0073C2
  "#EFC000FF", #EFC000
  "#ff9900"
)
names(protocol_color_sheme) <- protocols_int

save_ <- function(name, plot=NULL, data=NULL, width=NA, height=NA, dpi=600, use_pdf=TRUE, do_ggsave=TRUE){
  if (!is.null(plot)){
    plot %>% saveRDS(file.path(data_dir_, paste0(name, ".plot.rds")))
    format = ifelse(use_pdf, "pdf", "png")
    plot_path = file.path(figures_dir_, paste0(name, ".", format))
    print(plot_path)
    if (do_ggsave){
      ggsave(
        plot_path,
        plot,
        width=width,
        height=height,
        dpi=dpi
      )
    } else{
      do.call(format, list(
        plot_path,
        width=width,
        height=height
      ))
      print(plot)
      dev.off()
    }
  }
  if (!is.null(data)){
    data %>% saveRDS(file.path(data_dir_, paste0(name, ".rds")))
  }
}


read_ <- function(name, analysis_name, read_data=TRUE){
  if (read_data){
    print(file.path(data_dir, analysis_name, paste0(name, ".rds")))
    data <- readRDS(file.path(data_dir, analysis_name, paste0(name, ".rds"))) %>%
      return
  }
  else {
    print(file.path(data_dir, analysis_name, paste0(name, ".plot.rds")))
    readRDS(file.path(data_dir, analysis_name, paste0(name, ".plot.rds"))) %>% 
      return
  }
}

methrix_obj_dir = "/omics/groups/OE0219/internal/yuyu/benchmarking/analysis/read_in_hdf5_manuscript"

unified_pg <- theme(panel.background = element_rect(fill = "white", colour = "black", linetype="solid"), 
                    panel.border=element_rect(fill=NA), 
                    strip.background =element_blank())

mx_objs_name = c("GSNAP", "methylCtools", "bwameth", "gemBS", "BAT",
                 "Bismark", "FAME", "Biscuit", "BSBolt", "methylpy")

color_palette_color <- function(){
  return(scale_color_npg())
}

color_palette_fill <- function(){
  return(scale_fill_npg())
}


## RnBead cutoff
diffmeth_p_adj_fdr = 0.1

directions <- c("loss in tumor", "gain in tumor")  
