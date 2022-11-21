analysis_name <- "01-3_calc_depth_vs_coverage"
config = file.path(getwd(), "0_project_setting.R")
source(config)
library(methrix)

for(p in protocols_int){
  meth = methrix::load_HDF5_methrix(file.path(methrix_obj_dir, workflows_int[1]))
  meth = meth[,which(meth@colData$method==p)]

  for(c in workflows_int[2:length(workflows_int)]){
    meth_ = methrix::load_HDF5_methrix(file.path(methrix_obj_dir, c))
    meth_ = meth_[,which(meth_@colData$method==p)]
    meth = meth %>% combine_methrix(meth_, by="col")
  }

  cov_filt_tb <- meth %>%
    get_matrix(type="C", add_loci=T) %>%
    as_tibble

  save_(paste0("cov_filt_tb_",p), data=cov_filt_tb)
  
  rm(meth)

  beta_cov_tb <- cov_filt_tb %>%
    dplyr::rename(seqnames=chr) %>%
    dplyr::select(-strand) %>%
    gather(
      "sample",
      "cov",
      -seqnames, -start,
    ) %>%
    mutate(
      workflow=sapply(sample, function(s){
        str_split(s, "\\.")[[1]][4]
      })
    ) %>%
    mutate(
      sample=sapply(sample, function(s){
        str_split(s, "\\.")[[1]][3]
      })
    ) %>%
    select(sample, workflow, cov, everything()) %>%
    dplyr::rename(chr=seqnames) %>%
    dplyr::rename(pos=start)
  
  save_(paste0("beta_cov_tb_",p), data=beta_cov_tb)  
  rm(cov_filt_tb)
  rm(beta_cov_tb)
}
