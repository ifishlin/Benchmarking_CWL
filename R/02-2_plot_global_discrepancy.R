analysis_name <- "02-2_plot_global_discrepancy"
config = file.path(getwd(), "0_project_setting.R")
source(config)

m = read_("discrepancy_tb", "02-1_calc_global_discrepancy")

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

mat = mat[idx,idx]

colnames(mat)=gsub("^.*\\&","",colnames(mat))
rownames(mat)=gsub("^.*\\&","",rownames(mat))

library(corrplot)

######
library(viridis)
col <- viridis(200, direction =-1)
psize = ceiling(ncol(mat)/3)*100
pdf(height=16, width=16, pointsize=11, file=file.path(figures_dir_, "discrepancy_heatmap.pdf"))
corrplot.mixed(mat, is.corr=FALSE, 
               lower = 'number', 
               upper = 'color', 
               upper.col = col,
               lower.col = col,
               addCoefasPercent = TRUE,
               tl.col = "black",
               tl.pos = 'lt', 
               diag = 'u')
dev.off()


for(s in c("5N", "5T", "6N", "6T")){

  m = read_(paste0("discrepancy_",s,"_tb"), "02_Discrepancy_calculation")
  
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
  
  mat = mat[idx,idx]
  
  colnames(mat)=gsub("^.*\\&","",colnames(mat))
  rownames(mat)=gsub("^.*\\&","",rownames(mat))
  
  library(corrplot)
  
  ######
  library(viridis)
  col <- viridis(200, direction =-1)
  psize = ceiling(ncol(mat)/3)*100
  pdf(height=16, width=16, pointsize=11, file=file.path(figures_dir_, paste0("discrepancy_heatmap_",s ,".pdf")))
  corrplot.mixed(mat, is.corr=FALSE, 
                 lower = 'number', 
                 upper = 'color', 
                 upper.col = col,
                 lower.col = col,
                 addCoefasPercent = TRUE,
                 tl.col = "black",
                 tl.pos = 'lt', 
                 diag = 'u')
  dev.off()  
}
