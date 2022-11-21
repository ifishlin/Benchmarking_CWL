analysis_name <- "01-5_plot_beta_dist"
config = file.path(getwd(), "0_project_setting.R")
source(config)

## Load methrix
##########################
meth_list = lapply(workflows, function(x){
  methrix_o <- methrix::load_HDF5_methrix(file.path(methrix_obj_dir, x))
}) 

## Sampling
sample_size=500
set.seed(2021)
idx = sample(1:29401795, sample_size) # The number of CpGs is 29,401,795
combs=c()
df=data.frame(matrix(0, ncol = 0, nrow = sample_size))
for(obj in meth_list){
  combs = c(combs,paste(obj@colData@listData$sample, 
                        replace_prot_name_(obj@colData@listData$method),  
                        replace_wf_name_(obj@colData@listData$pipeline), sep = "_"))
  beta = obj@assays@data@listData[["beta"]][idx, ]
  df = cbind(df, as.data.frame(beta))
}
colnames(df) = combs
save_("beta_dist_sampling_tb", data=df)

## Reformat
df_tmp = df %>% 
  gather(key="MesureType", value="beta") %>%
  mutate(sample=gsub("(5N|5T|6T|6N)_(WGBS|T-WGBS|PBAT|Swift|EM-seq)_(.*)", "\\1", MesureType)) %>%
  mutate(protocol=gsub("(5N|5T|6T|6N)_(WGBS|T-WGBS|PBAT|Swift|EM-seq)_(.*)", "\\2", MesureType)) %>%
  mutate(workflow=gsub("(5N|5T|6T|6N)_(WGBS|T-WGBS|PBAT|Swift|EM-seq)_(.*)", "\\3", MesureType)) %>% 
  dplyr::select(-MesureType)
  
#df_tmp$workflow <- replace_wf_name_(df_tmp$workflow)
df_tmp$protocol =  factor(df_tmp$protocol, levels = protocols)
  
## Plot
p <- ggplot(df_tmp, aes(x=beta, color=workflow)) +
  geom_density() +   
  theme(panel.background = element_rect(fill = "white", colour = "black", linetype="solid"),
                           panel.border=element_rect(fill=NA),
                           strip.background = element_blank(),
                           legend.position='right') + 
  facet_wrap(~protocol, ncol=5, scales = "free_y") + color_palette_color() + color_palette_fill()

save_(
  "beta_dist_density",
  use_pdf=TRUE,
  plot=p,
  width=10,
  height=4
)

