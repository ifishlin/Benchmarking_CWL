analysis_name <- "02-4_plot_discrepancy_by_annot"
config = file.path(getwd(), "0_project_setting.R")
source(config)
library(annotatr)
library(methrix)
annots_hg38 = builtin_annotations()[str_detect(builtin_annotations(), "hg38")]
annots_hg38 = annots_hg38[-17]

b = lapply(annots_hg38, function(k){
  f = paste0(data_dir, "/02-3_calc_discrepancy_by_annot", paste0("/discrepancy_tb_", k, ".rds"))

  if(file.exists(f)){
    
    m = read_(paste0("discrepancy_tb_", k), "02-3_calc_discrepancy_by_annot")
    
    mat = as.matrix(m)
    ##### reorder by protocol
    idx = order(colnames(mat))
    mat = mat[idx,idx]
    
    colnames(mat) = replace_wf_name_(colnames(mat))
    rownames(mat) = replace_wf_name_(rownames(mat))
    
    ##### second reorder by workflow
    idx=c()
    idx=c(idx, which(str_detect(colnames(mat), "^WGBS"), TRUE))
    idx=c(idx, which(str_detect(colnames(mat), "SWIFT"), TRUE))
    idx=c(idx, which(str_detect(colnames(mat), "TWGBS"), TRUE))
    idx=c(idx, which(str_detect(colnames(mat), "PBAT"), TRUE))
    idx=c(idx, which(str_detect(colnames(mat), "EMSEQ"), TRUE))
    
    df_mat = data.frame(mat[idx,idx])
    
    mat = mat[idx,idx]
    
    #colnames(mat)=gsub("^.*\\&","",colnames(mat))
    #rownames(mat)=gsub("^.*\\&","",rownames(mat))
    
    mat_condensed = lapply(c("WGBS", "SWIFT", "TWGBS", "PBAT", "EMSEQ"), function(x){
      idx = which(str_detect(colnames(mat), paste0("^", x)), TRUE)
      a = df_mat[,idx]
      b = data.frame(rowMeans(a))  
    }) %>% bind_cols()   
    
    colnames(mat_condensed) = c("WGBS", "SWIFT", "TWGBS", "PBAT", "EMSEQ")
  
    mat_condensed = mat_condensed %>% dplyr::select("WGBS")
    
    mat_condensed = mat_condensed %>% mutate(name=rownames(mat_condensed))
    mat_condensed$name = gsub("^.*\\&","",rownames(mat_condensed))
    
    mat_condensed$from = gsub("\\&.*","",rownames(mat_condensed))
    
    a = mat_condensed %>% gather("to", "value", -name, -from) %>% mutate(annotation=k)
    
  }
}) %>% bind_rows()

b$name = factor(b$name, level=rev(workflows))
b$from <- replace_prot_name_(b$from)
b$from = factor(b$from, levels=protocols)

b$annotation = gsub("^hg38_", "", b$annotation)

b = b %>% mutate(s=as.character(round(value,3)))

## genome-wide
m = read_("discrepancy_tb", "02-2_plot_global_discrepancy")

mat = as.matrix(m)
##### reorder by protocol
idx = order(colnames(mat))
mat = mat[idx,idx]

colnames(mat) = replace_wf_name_(colnames(mat))
rownames(mat) = replace_wf_name_(rownames(mat))

##### second reorder by workflow
idx=c()
idx=c(idx, which(str_detect(colnames(mat), "^WGBS"), TRUE))
idx=c(idx, which(str_detect(colnames(mat), "SWIFT"), TRUE))
idx=c(idx, which(str_detect(colnames(mat), "TWGBS"), TRUE))
idx=c(idx, which(str_detect(colnames(mat), "PBAT"), TRUE))
idx=c(idx, which(str_detect(colnames(mat), "EMSEQ"), TRUE))

df_mat = data.frame(mat[idx,idx])

mat = mat[idx,idx]

#colnames(mat)=gsub("^.*\\&","",colnames(mat))
#ownames(mat)=gsub("^.*\\&","",rownames(mat))

mat_condensed = lapply(c("WGBS", "SWIFT", "TWGBS", "PBAT", "EMSEQ"), function(x){
  idx = which(str_detect(colnames(df_mat), paste0("^", x)), TRUE)
  a = df_mat[,idx]
  b = data.frame(rowMeans(a))  
}) %>% bind_cols()   

colnames(mat_condensed) = c("WGBS", "SWIFT", "TWGBS", "PBAT", "EMSEQ")

mat_condensed = mat_condensed %>% select("WGBS")

mat_condensed = mat_condensed %>% mutate(name=rownames(mat_condensed))
mat_condensed$name = gsub("^.*\\&","",rownames(mat_condensed))

mat_condensed$from = gsub("\\&.*","",rownames(mat_condensed))

h = mat_condensed %>% gather("to", "value", -name, -from) %>% mutate(annotation="whole_genome") %>% 
  mutate(s=as.character(round(value,3)))

h$name = factor(h$name, level=rev(workflows))
h$from <- replace_prot_name_(h$from)
h$from = factor(h$from, levels=protocols)

#
b = rbind(b, h)


ordered = b %>% filter(from=="WGBS") %>% group_by(annotation) %>% 
  summarise(m=mean(value)) %>% arrange(m) %>% select(annotation)

b$annotation = factor(b$annotation, levels=as.list(ordered)$annotation)

ordered_w = b %>% filter(from=="WGBS") %>% group_by(name) %>% 
  summarise(m=mean(value)) %>% arrange(m) %>% select(name)

b$name = factor(b$name, levels=rev(as.list(ordered_w)$name))

tmp = data.frame(name=c("BAT"), from="PBAT", to="WGBS", value=NA, annotation=as.list(ordered)$annotation, s=NA)
b2 = rbind(b, tmp)
tmp = data.frame(name=c("gemBS"), from="PBAT", to="WGBS", value=NA, annotation=as.list(ordered)$annotation, s=NA)
b2 = rbind(b2, tmp)

g = b2 %>% ggplot(aes(annotation, name, fill=value)) + 
  geom_tile(colour="white",size=0.25) +
  theme(
    panel.border=element_blank(),
    panel.background = element_blank(),
    legend.text=element_text(face="bold"),
    axis.ticks=element_line(size=0.4),
    axis.title.x = element_text(),
    strip.background =element_blank(),
    axis.text.x = element_text(angle = 60, hjust=1, size=7), 
    legend.position='bottom',
    strip.text = element_text(size = 11)
  ) + facet_wrap(~ from, ncol=5) + ylab("") + scale_fill_gradient(low = "yellow2",
                                                       high = "darkgreen",
                                                       guide = "colorbar")

g

save_("wg_annotation_plot", plot=g, use_pdf = TRUE, width=11, height=5)

