analysis_name <- "k04_meth_and_cov_distr"
source("/omics/groups/OE0219/internal/yuyu/10.Benchmarking/analysis/data_processing_redo/00.project_setting.R")
library("methrix")

# create methrix
# mcall_corr_mbias_dir="/omics/groups/OE0219/internal/yuyu/DNA_meth_protocol_benchmarking/processed/bismark/mbias_corr_mcall_resoftlink"
# bdg_files <- sapply(samples, function(s){
#   sapply(protocols_int, function(p){
#     file.path(mcall_corr_mbias_dir, paste0(s,"_", p, ".bismark.cov.gz"))
#   })
# })
# names(bdg_files) <- sapply(samples, function(s){
#   sapply(protocols_int, function(p){
#     paste0(s, "_", p)
#   })
# })
# 
# # filter non existing files:
# bdf_file_to_remove <- c()
# for (bdg_file in bdg_files){
#   if (! file.exists(bdg_file)){
#     bdf_file_to_remove <- c(bdf_file_to_remove, bdg_file)
#   }
# }
# bdg_files <- bdg_files[ ! bdg_files %in% bdf_file_to_remove ]
# 
# 
# sample_ids <- sapply(str_split(names(bdg_files), "_"), function(x) x[1])
# protocols_ <- sapply(str_split(names(bdg_files), "_"), function(x) x[2])
# unique_id <- paste0(protocols_, "_", sample_ids)
# 
# 
# sample_anno <- data.frame(
#   row.names = unique_id,
#   unique_id = unique_id,
#   sample = sample_ids,
#   protocol = protocols_,
#   tumor = substring(sample_ids, 2, 2)=="T",
#   patient = substring(sample_ids, 1, 1) %>% as.integer(),
#   assembly = rep("hg38", length(sample_ids)),
#   pipeline = rep("bismark", length(sample_ids)),
#   stringsAsFactors = FALSE
# )
# 
# ## methrix object:
# ref_cpgs = methrix::extract_CPGs(ref_genome = "BSgenome.Hsapiens.UCSC.hg38")
# 
# meth <- methrix::read_bedgraphs(
#   files = bdg_files,
#   pipeline="Bismark_cov",
#   zero_based=FALSE,
#   ref_cpgs = ref_cpgs,
#   coldata = sample_anno,
#   collapse_strands = TRUE,
#   stranded = T
# )
# 
# save_(
#   "meth",
#   data=meth
# )
# 
# unique_id_to_protocol <- sapply(unique_id, function(u){
#   strsplit(u, "_")[[1]][1]
# })
# 
# unique_id_to_sample <- sapply(unique_id, function(u){
#   strsplit(u, "_")[[1]][2]
# })
# 
# # as tibble coverage values:
# cov_tb <- meth %>%
#   get_matrix(type="C", add_loci=T) %>%
#   as_tibble %>%
#   select(-strand)
# 
# cov_tb <- cov_tb %>%
#   gather(
#     "sample",
#     "cov",
#     -chr,
#     -start
#   ) %>%
#   mutate(
#     protocol=unique_id_to_protocol[sample],
#     sample=unique_id_to_sample[sample]
#   ) %>%
#   dplyr::rename(
#     pos=start
#   ) %>%
#   select(sample, protocol, chr, pos, cov) %>%
#   mutate_at(vars(sample, protocol, chr), factor) %>%
#   mutate(pos=as.integer(pos), cov=as.integer(cov)) %>%
#   mutate(cov=ifelse(is.na(cov), 0, cov)) %>%
#   set_order()
# 
# save_(
#   "cov_tb",
#   data=cov_tb
# )

############# Functions for ploting ############# 

get_distribution_overview <- function(beta_cov_tb, var="cov", add_category=NULL, condense=FALSE){
  if (condense){
    dist_overview_tb <- beta_cov_tb %>%
      group_by(protocol) 
  }
  else {
    dist_overview_tb <- beta_cov_tb %>%
      group_by(sample, protocol) 
  }
  dist_overview_tb <- dist_overview_tb %>%
    dplyr::summarize(
      mean=mean(UQ(sym(var)), na.rm=T),
      sd=sd(UQ(sym(var)), na.rm=T),
      med=median(UQ(sym(var)), na.rm=T),
      mad=mad(UQ(sym(var)), na.rm=T),
      upper=quantile(UQ(sym(var)), 0.75, na.rm=T),
      lower=quantile(UQ(sym(var)), 0.25, na.rm=T),
      n=n()
    ) %>%
    ungroup() %>%
    dplyr::mutate(se=sd/n)
  if (!is.null(add_category)){
    dist_overview_tb <- dist_overview_tb %>%
      mutate(category=add_category) %>%
      dplyr::select(category, everything())
  }
  return(dist_overview_tb)
}

plot_1_minus_ecdf_cov <- function(
    beta_cov_tb, add_category=NULL, condense=FALSE, facet=TRUE, x_limits=c(1, 100), log_scale=TRUE
){
  protocol_ypos <- seq(-0.25, -0.05, length.out=5)
  names(protocol_ypos) <- protocols
  med_cov_tb <-  get_distribution_overview(beta_cov_tb, "cov", add_category, condense=condense) %>%
    dplyr::mutate(ypos=protocol_ypos[protocol])
  
  save_(paste0("med_cov_tb_",condense),data=med_cov_tb)
  
  if(facet || condense){
    plot_ <- beta_cov_tb %>%
      ggplot() +
      geom_line(
        aes(x=cov, y=1 - ..y.., color=protocol), 
        stat="ecdf",
        alpha=0.7,
        size=1
      ) +
      geom_point(
        aes(med, ypos, color=protocol),
        size=3.5,
        alpha=1,
        shape=15,
        data=med_cov_tb
      ) +
      geom_segment(
        aes(
          x=lower, xend=upper, y=ypos, yend=ypos,
          color=protocol
        ),
        size=3.5,
        alpha=0.5,
        data=med_cov_tb
      ) +
      scale_y_continuous(
        breaks=c(-0.15, seq(0,1,0.25)), 
        labels=c("median\n & IQR", 0.00, 0.25, 0.5, 0.75, 1.00),
        limits=c(-0.25,1)
      )
  }
  else {
    plot_ <- beta_cov_tb %>%
      ggplot() +
      geom_line(
        aes(x=cov, y=1 - ..y.., color=protocol, alpha=sample, linetype=sample), 
        stat="ecdf",
        alpha=0.7,
        size=1
      ) +
      scale_alpha_manual(values=sample_alpha_sheme) +
      scale_linetype_manual(values=sample_linetype_sheme)
  }
  plot_ <- plot_ +
    ylab("1 - ecdf") +
    xlab("coverage") +
    scale_color_manual(values=protocol_color_sheme)
  
  if(log_scale){
    plot_ <- plot_ +
      scale_x_log10(breaks=c(1, 10, 15, 30, 50, 100), limits=x_limits)
  }
  
  if(facet){
    plot_ <- plot_ + 
      facet_grid(. ~ sample) +
      theme(
        panel.border=element_rect(fill=NA),
        strip.background =element_blank()
        #theme(legend.position = "none")
      )
  }
  plot_ <- plot_ +
  theme(panel.background = element_rect(fill = "white", colour = "black", linetype="solid"),
        panel.border=element_rect(fill=NA),
        strip.background = element_blank(),
        legend.position='none')      
  
  return(
    list(
      data=med_cov_tb,
      plot=plot_
    )
  )
}

beta_cov_tb <- read_(
  "cov_tb",
  "k04_meth_and_cov_distr"
) 

replace_prot_name_ <- function(protocol_name){
  protocol_name <- gsub("SWIFT", "Swift", protocol_name)
  protocol_name <- gsub("TWGBS", "T-WGBS", protocol_name)
  return(protocol_name)
}

protocol_color_sheme <- c(
  "#CD534CFF",
  "#33cc00",
  "#0073C2FF",
  "#EFC000FF",
  "#ff9900"
)

protocols_show <- c(
  "WGBS",
  "Swift",
  "T-WGBS",
  "PBAT",
  "EMseq"
)

names(protocol_color_sheme) <- protocols_show

beta_cov_tb$protocol <- replace_prot_name_(beta_cov_tb$protocol)
beta_cov_tb$protocol = factor(beta_cov_tb$protocol, levels = protocols_show)

one_minus_ecdf_cov_condensed_plot <- plot_1_minus_ecdf_cov(beta_cov_tb)$plot
#
# ## what I need
save_(
  "1_minus_ecdf_coverage",
  plot=one_minus_ecdf_cov_condensed_plot,
  width=7.5,
  height=3.5
)

save_(
  "1_minus_ecdf_coverage",
  plot=one_minus_ecdf_cov_condensed_plot,
  use_pdf = TRUE,
  width=7.5,
  height=3.5
)

one_minus_ecdf_cov_plot <- plot_1_minus_ecdf_cov(beta_cov_tb, condense=TRUE, facet=FALSE)$plot
# ## what I need
save_(
  "1_minus_ecdf_coverage_condense",
  plot=one_minus_ecdf_cov_plot,
  width=6,
  height=5
)

save_(
  "1_minus_ecdf_coverage_condense",
  plot=one_minus_ecdf_cov_plot,
  use_pdf = TRUE,
  width=6,
  height=5
)