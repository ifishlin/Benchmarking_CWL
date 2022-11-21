analysis_name <- "06-1_plot_ranking_heatmap"
config = file.path(getwd(), "0_project_setting.R")
source(config)

##### idx
wgbs_idx  =  1:10
swift_idx = 11:20
twgbs_idx = 21:30
pbat_idx  = 31:38
emseq_idx = 39:48 
# idx for metrics where some workflows failed in PBAT. 
idxs = list(wgbs_idx, swift_idx, twgbs_idx, pbat_idx, emseq_idx)
# idx for metrics 
idxs_all = list(wgbs_idx, swift_idx, twgbs_idx, c(31:40), c(41:50))

##### 1. whole genome 
m = read_("discrepancy_tb", "02-1_calc_global_discrepancy")

##### data transformation
mat = as.matrix(m)

##### first reorder by protocol
idx = order(colnames(mat))
mat = mat[idx,idx]

##### second reorder by workflow
idx=c()
idx=c(idx, which(str_detect(colnames(mat), "^WGBS"), TRUE))
idx=c(idx, which(str_detect(colnames(mat), "SWIFT|Swift"), TRUE))
idx=c(idx, which(str_detect(colnames(mat), "TWGBS|T-WGBS"), TRUE))
idx=c(idx, which(str_detect(colnames(mat), "PBAT"), TRUE))
idx=c(idx, which(str_detect(colnames(mat), "EMSEQ|EM-seq"), TRUE))
df_mat = data.frame(mat[idx,idx])

##### compare to the WGBS profiles
df_mat = df_mat[,1:10]
df_rsums = data.frame(rowSums(df_mat))
df_rsums$workflow = gsub("(.*)&(.*)", "\\2", rownames(df_rsums))
df_rsums$protocol = gsub("(.*)&(.*)", "\\1", rownames(df_rsums))
colnames(df_rsums) <- c("score", "workflow", "protocol")

df_rsums$protocol <- replace_prot_name_(df_rsums$protocol)
df_rsums$workflow <- replace_wf_name_(df_rsums$workflow)
rownames(df_rsums) = NULL

#####
tb1 = data.frame(score=c(), workflow=c(), protocol=c())
for(idx in idxs){

  m = df_rsums[idx,] 
  mean = mean(m$score)
  sd = sd(m$score)
  maxv = max(m$score)
  minv = min(m$score)
  d = maxv - minv
  
  m <- m %>% mutate(zscore=(score-mean)/sd, 
                    minmax=(score-minv)/d, 
                    mtype="Profile SI")

  s = sort(m$score, index=TRUE)$ix
  m = cbind(m, rank=s, annot = s)
  
  tb1 = rbind(tb1, m)
}

## 2. gold-standard
beta_cov_gs_ds_tb <- read_("beta_cov_gs_tb", "03-1_calc_beta_cov_gs")
beta_cov_gs_ds_tb$protocol <- replace_prot_name_(beta_cov_gs_ds_tb$protocol)
beta_cov_gs_ds_tb$protocol <- factor(beta_cov_gs_ds_tb$protocol, levels=protocols)
beta_cov_gs_ds_tb$workflow <- replace_wf_name_(beta_cov_gs_ds_tb$workflow)

beta_cov_gs_ds_tb <- beta_cov_gs_ds_tb %>%
  mutate(dev=apply(cbind(beta, lower, upper), 1, 
                   function(r){
                     m=r[1]
                     lower=r[2]
                     upper=r[3]
                     return(
                       ifelse(m>=lower, ifelse(m<=upper, 0, m-upper), m-lower) 
                     )
                   })
  )

m2 = beta_cov_gs_ds_tb %>% dplyr::select("protocol", "workflow", "dev")
m2 = m2[!is.na(m2$dev),]
m2 = m2 %>% dplyr::group_by(protocol, workflow) 
m2 = m2 %>% dplyr::summarise(score = mean(abs(dev))) %>% dplyr::arrange(protocol)

tb2 = data.frame()
for(idx in idxs){

  m = m2[idx,] 
  
  mean = mean(m$score)
  sd = sd(m$score)
  maxv = max(m$score)
  minv = min(m$score)
  d = maxv - minv
  
  m <- m %>% mutate(zscore=(score-mean)/sd, 
                    minmax=(score-minv)/d,
                    mtype="GS dev")

  s = sort(m$score, index=TRUE)$ix
  m = cbind(m, rank=s, annot = s)
  
  tb2 = rbind(tb2, m)
}

## 3. DMI AUC
ROCit_tb = read_("weighted_AUC", "04-7_plot_weighted_AUC")
df_auc = unique(ROCit_tb %>% dplyr::select("workflow", "protocol", "weighted"))
df_auc <- df_auc %>% dplyr::arrange(protocol) %>% dplyr::rename(score = weighted)
rownames(df_auc) <- NULL

df_auc$protocol <- replace_prot_name_(df_auc$protocol)
df_auc$workflow <- replace_wf_name_(df_auc$workflow)
df_auc$protocol <- factor(df_auc$protocol, levels=protocols)
df_auc <- df_auc %>% dplyr::arrange(protocol)

tb3 = data.frame()
for(idx in idxs){

  m = df_auc[idx,] 

  mean = mean(m$score)
  sd = sd(m$score)
  maxv = max(m$score)
  minv = min(m$score)
  d = maxv - minv
  
  m <- m %>% mutate(zscore=(score-mean)/sd, 
                    minmax=(score-minv)/d,
                    mtype="DMI AUC")
 
  s = sort(m$score, index=TRUE, decreasing=TRUE)$ix
  m = cbind(m, rank=s, annot = s)
  
  tb3 = rbind(tb3, m)
}

## 4. WG_cov_ncpg
###############

tb = read_("wg_cov_ncpg_condesed_tb", "01-2_plot_CpGscovered_depth") %>% arrange(protocols)
tb = tb %>% ungroup %>% dplyr::mutate(score=ifelse(s<0, -s, s), workflow=workflows, protocol=protocols) %>% 
  dplyr::select(-c(s, workflows, protocols))
tb$protocol <- replace_prot_name_(tb$protocol)
tb_ncpg = tb %>% dplyr::filter(mtype=="n_cpg") %>% dplyr::arrange(protocol)
tb_cov  = tb %>% dplyr::filter(mtype=="coverage") %>% dplyr::arrange(protocol)

tb4 = data.frame()
for(idx in idxs){

  m = tb_ncpg[idx,]

  mean = mean(m$score)
  sd = sd(m$score)
  maxv = max(m$score)
  minv = min(m$score)
  d = maxv - minv
  
  m <- m %>% mutate(zscore=(score-mean)/sd, 
                    minmax=(score-minv)/d,
                    mtype="n_cpg")
  
  s = rev(sort(m$score, index=TRUE)$ix)
  m = cbind(m, rank=s, annot = s)
  
  tb4 = rbind(tb4, m)  
}

tb5 = data.frame()
for(idx in idxs){
  m = tb_cov[idx,]
  
  mean = mean(m$score)
  sd = sd(m$score)
  maxv = max(m$score)
  minv = min(m$score)
  d = maxv - minv
  
  m <- m %>% mutate(zscore=(score-mean)/sd,
                    minmax=(score-minv)/d,
                    mtype="coverage")

  s = rev(sort(m$score, index=TRUE)$ix)

  m = cbind(m, rank=s, annot = s) 
  
  tb5 = rbind(tb5, m)  
}

## hw_performance
###############
hw_tb = read_("hw_runt_mem_na_tb", "05-1_plot_runtime_and_maxmem")
hw_tb$protocol <- replace_prot_name_(hw_tb$protocol)
hw_tb$protocol <- factor(hw_tb$protocol, levels=protocols)
hw_tb$workflow <- replace_wf_name_(hw_tb$workflow)
hw_run_time_tb <- hw_tb %>% dplyr::mutate(score=run_time_h) %>% dplyr::arrange(protocol) %>% dplyr::select(protocol, workflow, score)
hw_max_mem_tb <- hw_tb %>% dplyr::mutate(score=max_mem_g) %>% dplyr::arrange(protocol) %>% dplyr::select(protocol, workflow, score)

##

tb6 = data.frame()
for(idx in idxs){

  m = hw_run_time_tb[idx,]
  
  mean = mean(m$score)
  sd = sd(m$score)
  maxv = max(m$score)
  minv = min(m$score)
  d = maxv - minv
  
  m <- m %>% mutate(zscore=(score-mean)/sd,
                    minmax=(score-minv)/d,
                    mtype="run_t")
  
  s = rev(sort(m$score, index=TRUE)$ix)
  
  m = cbind(m, rank=s, annot = s) 
  
  tb6 = rbind(tb6, m)  
}

tb7 = data.frame()
for(idx in idxs){
  m = hw_max_mem_tb[idx,]
  
  mean = mean(m$score)
  sd = sd(m$score)
  maxv = max(m$score)
  minv = min(m$score)
  d = maxv - minv
  
  m <- m %>% mutate(zscore=(score-mean)/sd,
                    minmax=(score-minv)/d,
                    mtype="max_mem")
  
  s = rev(sort(m$score, index=TRUE)$ix)
  
  m = cbind(m, rank=s, annot = s) 
  
  tb7 = rbind(tb7, m) 
}
##
###############
tb_ranking = rbind(tb1, tb2, tb3, tb4, tb5, tb6, tb7)
tb_ranking$mtype = factor(tb_ranking$mtype, levels=c('coverage', 'n_cpg', "Profile SI", "GS dev", "DMI AUC", "run_t", "max_mem"))

save_("final_ranking_normalization_tb", data=tb_ranking)

# excel output
tb_exl = data.frame()
for(pro in protocols){
  a = tb_ranking %>% filter(protocol==pro) %>% dplyr::select(workflow, mtype, score)
  a = reshape(a, idvar="workflow", timevar="mtype", direction="wide") %>% mutate(protocol=pro)
  tb_exl = rbind(tb_exl, a)
}
tb_exl = tb_exl %>% dplyr::select(workflow, protocol, everything())

tb_exl %>%
  writexl::write_xlsx(path = file.path(data_dir_, "score_table_run.xls"))
#

tb_ranking$protocol = factor(tb_ranking$protocol, levels=protocols)

## create blank rows for unavailable 
tb_lack = data.frame(score=NA, workflow=c('gemBS'), protocol=c('PBAT'), 
                     rank=(length(pbat_idx)+1)/2, zscore=NA, minmax=NA, annot=NA,
                     mtype=c('coverage', 'n_cpg', 'Profile SI', 'GS dev', 'DMI AUC', 'run_t', "max_mem"))
tb_ranking = rbind(tb_ranking, tb_lack)
tb_lack = data.frame(score=NA, workflow=c('BAT'), protocol=c('PBAT'), 
                     rank=(length(pbat_idx)+1)/2, zscore=NA, minmax=NA, annot=NA,
                     mtype=c('coverage', 'n_cpg', 'Profile SI', 'GS dev', 'DMI AUC', 'run_t', "max_mem"))
tb_ranking = rbind(tb_ranking, tb_lack)

save_("final_ranking_normalization_tb", data=tb_ranking)

## order by mean rank
global_ranking = tb_ranking %>% dplyr::group_by(workflow) %>% 
    dplyr::summarize(s_=mean(as.numeric(rank))) %>%
    dplyr::arrange(s_, workflow) 

## order workflow by global_ranking
tb_ranking$workflow = factor(tb_ranking$workflow, levels=rev(global_ranking$workflow))

## Label only top three
tb_ranking_ = tb_ranking %>% dplyr::mutate(sn_ = ifelse(annot>3, NA, annot)) %>% 
  dplyr::mutate(annot=as.character(annot))

tb_ranking_$annot = factor(tb_ranking_$annot, levels=1:11)

save_("final_ranking_normalization_data_tb", data=tb_ranking_)

mPalette <- c("#004616", "#005A1E", "#007328",
              "#007F2E", "#2A964D", "#59A76C",
              "#B2D3B9", "#CCE1D0", "#DFEBE2", "#ECF2ED", "darkgray")

g = tb_ranking_ %>% 
  ggplot(aes(mtype, workflow, fill= annot)) + 
  geom_tile(colour="white",size=0.25) + geom_text(aes(label = sn_), color="white") +
  facet_wrap(~ protocol, ncol=5) +
  theme(
    panel.border=element_blank(),
    panel.background = element_blank(),
    legend.text=element_text(face="bold"),
    axis.ticks=element_line(size=0.4),
    axis.title.x = element_text(),
    strip.background =element_blank(),
    axis.text.x = element_text(angle = 60, hjust=1), legend.position='bottom'
  ) + guides(fill=guide_legend(title="Rank",ncol = 11, byrow = TRUE)) + xlab("") +  ylab("") + coord_fixed() +
  scale_x_discrete(breaks=c('coverage', 'n_cpg', 'Profile SI', 'GS dev', 'DMI AUC', "run_t", "max_mem"), 
                   label = c("Depth", '% CpG covered', 'Profile_SIM', 'GS_DEV', 'DML_wAUC', "Run_Time", "Max_Mem")) + 
  scale_colour_manual(values=mPalette) + 
  scale_fill_manual(values=mPalette)

g

save_(paste0("heatmap_grid_ten_gsmean_rank"), plot=g, use_pdf = TRUE, width=8, height=5)

g = tb_ranking_ %>% 
  ggplot(aes(mtype, workflow, fill= annot)) + 
  geom_tile(colour="white",size=0.25) + geom_text(aes(label = annot), color="white") +
  facet_wrap(~ protocol, ncol=5) +
  theme(
    panel.border=element_blank(),
    panel.background = element_blank(),
    legend.text=element_text(face="bold"),
    axis.ticks=element_line(size=0.4),
    axis.title.x = element_text(),
    strip.background =element_blank(),
    axis.text.x = element_text(angle = 60, hjust=1), legend.position='bottom'
  ) + guides(fill=guide_legend(title="Rank",ncol = 11, byrow = TRUE)) + xlab("") +  ylab("") + coord_fixed() +
  scale_x_discrete(breaks=c('coverage', 'n_cpg', 'Profile SI', 'GS dev', 'DMI AUC', "run_t", "max_mem"), 
                   label = c("Depth", '% CpG covered', 'Profile_SIM', 'GS_DEV', 'DML_wAUC', "Run_Time", "Max_Mem")) + 
  scale_colour_manual(values=mPalette) + 
  scale_fill_manual(values=mPalette)

g

save_(paste0("heatmap_grid_ten_gsmean_rank_10"), plot=g, use_pdf = TRUE, width=8, height=5)


## alluvial plot
indv_ranking_tb = data.frame(1:10)

p = "WGBS"
indv_ranking = tb_ranking %>% dplyr::filter(protocol==p) %>% group_by(workflow) %>%
  dplyr::summarize(s_=mean(as.numeric(rank))) %>%
  dplyr::arrange(s_) %>% ungroup %>% 
  mutate(r1=as.numeric(as.factor(rank(s_)))) 

indv_rank = as.character(indv_ranking$workflow)

indv_ranking = indv_ranking %>% 
               arrange(desc(workflow)) %>% 
               mutate(r1_n = indv_rank)
  
indv_ranking_tb = cbind(indv_ranking_tb, indv_ranking)

for(p in c("Swift", "T-WGBS", "PBAT", "EM-seq")){
  indv_ranking = tb_ranking %>% dplyr::filter(protocol==p) %>% group_by(workflow) %>%
    dplyr::summarize(s_=mean(as.numeric(rank))) %>%
    dplyr::arrange(s_) %>% ungroup %>% 
    mutate(r1=as.numeric(as.factor(rank(s_)))) 
  
  indv_rank = as.character(indv_ranking$workflow)
  
  indv_ranking = indv_ranking %>% 
    arrange(desc(workflow)) %>% mutate(r1_n = indv_rank) %>%
    dplyr::select(s_, r1, r1_n)

  indv_ranking_tb = cbind(indv_ranking_tb, indv_ranking)
}
indv_ranking_tb = indv_ranking_tb[-1]
colnames(indv_ranking_tb) = c("workflow", "WGBS_score", "WGBS_rank", "WGBS_rank_n", 
                              "Swift_score", "Swift_rank", "Swift_rank_n", 
                              "T-WGBS_score", "T-WGBS_rank", "T-WGBS_rank_n",  
                              "PBAT_score", "PBAT_rank", "PBAT_rank_n",  
                              "EM-seq_score", "EM-seq_rank", "EM-seq_rank_n") 

save_("indv_ranking_tb", data=indv_ranking_tb)


