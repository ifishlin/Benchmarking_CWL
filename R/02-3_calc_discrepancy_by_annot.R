analysis_name <- "02-3_calc_discrepancy_by_annot"
config = file.path(getwd(), "0_project_setting.R")
source(config)

library(methrix)
library(annotatr)

meth_list = lapply(workflows_int, function(x){
  methrix_o <- methrix::load_HDF5_methrix(file.path(methrix_obj_dir, x))
}) 

# workflows_int_ = sort(workflows_int)
# m_obj <- methrix::load_HDF5_methrix(paste0(methrix_obj_dir,"/",workflows_int_[1]))
# obj_list = list(m_obj)
# for(w in workflows_int_[2:length(workflows_int_)]){
#   m_obj <- methrix::load_HDF5_methrix(paste0(methrix_obj_dir,"/",w))
#   obj_list = append(obj_list, m_obj)
# }

# build matrix layout
pairs = c()
for(o in meth_list){
  o = o[1]
  combs = unique(paste(o@colData$method, o@colData$pipeline, sep="&"))
  pairs = c(pairs, combs)
}

pairs = sort(pairs)

# Select annotations for intersection with regions
# Note inclusion of custom annotation, and use of shortcuts
annots_hg38 = builtin_annotations()[str_detect(builtin_annotations(), "hg38")]
annots_hg38 = annots_hg38[-17] # remove hg38_lncrna_gencode

legal_chrs = c("chr1","chr2","chr3","chr4", "chr5", "chr6", "chr7","chr8","chr9","chr10",
              "chr11","chr12","chr13","chr14", "chr15", "chr16", "chr17","chr18","chr19","chr20",
              "chr21","chr22","chrX","chrY")

# List of annotation
# [1] "hg38_genes_1to5kb"               "hg38_genes_promoters"            "hg38_genes_cds"                  "hg38_genes_5UTRs"               
# [5] "hg38_genes_exons"                "hg38_genes_firstexons"           "hg38_genes_introns"              "hg38_genes_intronexonboundaries"
# [9] "hg38_genes_exonintronboundaries" "hg38_genes_3UTRs"                "hg38_genes_intergenic"           "hg38_cpg_islands"               
# [13] "hg38_cpg_shores"                 "hg38_cpg_shelves"                "hg38_cpg_inter"                  "hg38_enhancers_fantom"          
# [17] "hg38_basicgenes"                 "hg38_cpgs"  

for(k in annots_hg38){
  outfile_name = paste0(data_dir_,"/",paste0("discrepancy_tb","_", k, ".rds"))
  if(!file.exists(outfile_name)){
    next
  }
  
  print(k)
  annotations = build_annotations(genome = 'hg38', annotations = k)
  annotations = data.frame(annotations) %>% filter(seqnames %in% legal_chrs) %>% makeGRangesFromDataFrame(keep.extra.columns = T)

  df = data.frame(matrix(0, length(pairs), length(pairs)))
  rownames(df) = pairs
  colnames(df) = pairs
  
  df_5N = data.frame(matrix(0, length(pairs), length(pairs)))
  df_5T = data.frame(matrix(0, length(pairs), length(pairs)))
  df_6N = data.frame(matrix(0, length(pairs), length(pairs)))
  df_6T = data.frame(matrix(0, length(pairs), length(pairs)))
  rownames(df_5N) = colnames(df_5N) = pairs
  rownames(df_5T) = colnames(df_5T) = pairs
  rownames(df_6N) = colnames(df_6N) = pairs
  rownames(df_6T) = colnames(df_6T) = pairs
  
  df_samples = data.frame()  
  
  for(i in pairs[1:length(pairs)]){
    for(j in pairs[1:length(pairs)]){
      pair1 <- strsplit(i, split="&")
      p1 = pair1[[1]][1]
      w1 = pair1[[1]][2]
      idx_1 = match(w1, workflows_int)
      df_idx_1 = match(i, pairs)
      
      pair2 <- strsplit(j, split="&")
      p2 = pair2[[1]][1]
      w2 = pair2[[1]][2]
      idx_2 = match(w2, workflows_int)
      df_idx_2 = match(j, pairs)
      
      if(df_idx_1 >= df_idx_2) next
      
      obj_1 = meth_list[[idx_1]]
      m_annote = methrix::subset_methrix(m = obj_1[,which(obj_1@colData$method==p1)], regions=annotations)
      pair1_mat = methrix::get_matrix(m = m_annote)
      
      obj_2 = meth_list[[idx_2]]   
      m_annote = methrix::subset_methrix(m = obj_2[,which(obj_2@colData$method==p2)], regions=annotations)
      pair2_mat = methrix::get_matrix(m = m_annote)    
      difference <- pair1_mat-pair2_mat
      if (class(difference)=="DelayedMatrix"){
        sum_diff <- DelayedMatrixStats::colMeans2(abs(difference), na.rm = T)
      }else{
        sum_diff <- matrixStats::colMeans2(abs(difference), na.rm = T)
      }    

      df_5N[df_idx_1, df_idx_2] = sum_diff[1]
      df_5N[df_idx_2, df_idx_1] = sum_diff[1]
      df_5T[df_idx_1, df_idx_2] = sum_diff[2]
      df_5T[df_idx_2, df_idx_1] = sum_diff[2]
      df_6N[df_idx_1, df_idx_2] = sum_diff[3]
      df_6N[df_idx_2, df_idx_1] = sum_diff[3]
      df_6T[df_idx_1, df_idx_2] = sum_diff[4]
      df_6T[df_idx_2, df_idx_1] = sum_diff[4]
      
      df[df_idx_1, df_idx_2] = sum(sum_diff)/4
      df[df_idx_2, df_idx_1] = sum(sum_diff)/4
      print(paste0("++", i,"==",j," - ", sum(sum_diff)/4, "/", df_idx_1, "/", df_idx_2))
      
      df_samples = rbind(df_samples, c(i,j,sum_diff))
    }
  }
  
  colnames(df) = rownames(df) = replace_wf_name_(colnames(df))
  colnames(df_5N) = rownames(df_5N) = replace_wf_name_(colnames(df_5N))
  colnames(df_5T) = rownames(df_5T) = replace_wf_name_(colnames(df_5T))
  colnames(df_6N) = rownames(df_6N) = replace_wf_name_(colnames(df_6N))
  colnames(df_6T) = rownames(df_6T) = replace_wf_name_(colnames(df_6T))
    
  colnames(df_samples) = c("c1", "c2", "D5N", "D5T", "D6N","D6T")
  df_samples$D5N = as.numeric(df_samples$D5N)
  df_samples$D5T = as.numeric(df_samples$D5T)
  df_samples$D6N = as.numeric(df_samples$D6N)
  df_samples$D6T = as.numeric(df_samples$D6T)
  
  save_(paste0("discrepancy_tb","_",k), data=df)
  save_(paste0("discrepancy_sample_tb","_",k), data=df_samples)
  save_(paste0("discrepancy_5N_tb","_",k), data=df_5N)
  save_(paste0("discrepancy_5T_tb","_",k), data=df_5T)
  save_(paste0("discrepancy_6N_tb","_",k), data=df_6N)
  save_(paste0("discrepancy_6T_tb","_",k), data=df_6T)
}



