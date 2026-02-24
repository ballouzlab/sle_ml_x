library(dplyr)
library(tidyr)
library(tibble)
library(stringr)
library(ggplot2)
library(reshape2)
library(UpSetR)
library(ComplexHeatmap)
library(circlize)
library(ggsignif)
library(gridExtra)
library(ggpubr)
library(speckle)  
library(forcats)
library(beeswarm)

gene.set.colours <- c('chrX'='#8A0798', 'SLE'='#D90750', 'autosome' ='#D99750', 
                      'HVG.autosome' ='#D99799',  'HVG' ='#399750' )
method.colours <- c('boruta'='#BFFB00', 'enet'='#B875B1', 'intersection'='#D2EDF6', 'combined'='#4DB748')
model.colours <- c('logit'='#F4CE03', 'RF'='#BCEA9D', 'SVM'='#99B2F5', 'GBM'='#F5B29E', 'MLP'='#26779E', 'ensemble'='#F5A2F5')
gene.set = names(gene.set.colours )

#source('functions/replace.names.R')

#source('replace.names.R')
load("metadata_pbmc_female.control_managed.Rdata")
X.immune <- read.delim('X.immune.txt', header=FALSE)$V1
load('escapees.Rdata')
load('celltype.colours.RData')
celltype.colours <- colours

#katsir <- read.delim('Katsir.escape.txt')$Gene.Symbol
SLE <- read.delim('SLE_DisGeNet.tsv')$Gene
chrX <- read.delim('chrX_biomaRt.txt')
chrX <- subset(chrX, Gene.name != '')
chrX <- chrX$Gene.name



methods <- c('boruta', 'enet', 'intersection', 'combined')
models <- c('logit', 'RF', 'SVM', 'GBM', 'MLP')


# Read in metrics for all ML models across gene sets and splits
models_metrics_list <- list()
for(i in 1:10){
    metric.files <- unlist(lapply(methods, function(method){
    list.files(paste0('perez_update/split_', i, '/', method, '/metrics'), pattern='metrics_', full.names=TRUE)}))
    #list.files(paste0('pseudobulk/split_', i, '/', method, '/metrics'), pattern='metrics_', full.names=TRUE)}))
    # Check if metric files are empty
    metric.files <- metric.files[file.size(metric.files) > 0]
    metrics <- lapply(metric.files, read.csv)
    names(metrics) <- gsub('split_\\d+/|metrics/|.csv', '', metric.files) %>% 
        gsub('_metrics', '', .) %>%
        gsub('/', '_', .)
    metrics_df <- bind_rows(metrics, .id = 'celltype') %>%
    mutate(
        model = str_extract(celltype, "logit|RF|SVM|GBM|MLP"),
        gene.set = str_extract(celltype, "HVG\\.autosome|chrX|autosome|HVG|SLE"),
        method = str_extract(celltype, "boruta|enet|intersection|combined"),
        celltype = gsub("^perez_update_", "", celltype),
        #celltype = gsub("^pseudobulk_", "", celltype),
        celltype = gsub("^(boruta_|enet_|intersection_|combined_)?(logit_|RF_|SVM_|GBM_|MLP_)", "", celltype),
        celltype = gsub("\\.HVG\\.autosome|\\.SLE|\\.chrX|\\.autosome|\\.HVG", "", celltype),
        celltype = gsub("_", " ", celltype)
    )
    # Filter out autosome, HVG, and HVG.autosome gene sets
    #metrics_df <- subset(metrics_df, gene.set %in% c('chrX', 'SLE'))
    models_metrics_list[[i]] <- metrics_df
}

celltypes = unique(metrics_df$celltype)
temp = read.table("celltypes_test.txt", header=T, sep="\t") 
o = order(temp[,2]  ) 
temp = temp[o,]
celltypes.colours = turbo(23)
temp[,3] = celltypes.colours  
names(celltypes.colours) = celltypes[o] 
celltypes = celltypes[o] 

# Combine metrics from all splits
models_metrics_df <- bind_rows(models_metrics_list, .id='split')
models_metrics_df$model <- factor(models_metrics_df$model, levels=c('logit', 'RF', 'SVM', 'GBM', 'MLP'))
models_metrics_df$method <- factor(models_metrics_df$method, levels=c('boruta', 'enet', 'intersection', 'combined'))
models_metrics_df$gene.set <- factor(models_metrics_df$gene.set)
models_metrics_df$celltype <- factor(models_metrics_df$celltype, levels = celltypes)
models_metrics_df$split <- factor(models_metrics_df$split, levels=paste('split', 1:10, sep='_'))

#save(models_metrics_df, file='models_metrics_df.RData')

save(models_metrics_df, file='models_metrics_df_update.RData')

average_model_metrics_list = list() 
for (methodi in methods){  
average_model_metrics_df <- models_metrics_df %>% 
  subset(method == methodi)  %>%
  group_by(celltype, gene.set, model) %>%
  summarise(
    Accuracy = mean(Accuracy),
    Precision = mean(Precision),
    Recall= mean(Recall),
    F1=mean(F1),
    F1_lower=mean(F1_lower),
    F1_upper=mean(F1_upper),
    AUC=mean(AUC),
    AUC_lower=mean(AUC_lower),
    AUC_upper=mean(AUC_upper),
    AUPRC=mean(AUPRC),
    AUPRC_lower=mean(AUPRC_lower),
    AUPRC_upper=mean(AUPRC_upper),
    Kappa=mean(Kappa),
    MCC=mean(MCC),
    MCC_lower=mean(MCC_lower),
    MCC_upper=mean(MCC_upper),
    n_features=round(mean(n_features),0)) %>%
  data.frame()

  average_model_metrics_df$model <- factor(average_model_metrics_df$model, levels=c('logit', 'RF', 'SVM', 'GBM', 'MLP'))
  average_model_metrics_df$gene.set <- factor(average_model_metrics_df$gene.set)
  average_model_metrics_list[[methodi]] = average_model_metrics_df
  
} 
 

average_model_metrics_df = average_model_metrics_list$combined  
comparisons <- list(c('logit', 'MLP'), c('RF', 'MLP'), c('SVM', 'MLP'), c('GBM', 'MLP'))
ggplot(average_model_metrics_df  , aes(x=model, y=MCC, colour=model)) +
  geom_jitter(width=0.2) +
  geom_boxplot(outlier.shape = NA, color = 'black', fill = NA) +
  geom_signif(comparisons=comparisons, map_signif_level=TRUE, 
              test='wilcox.test', color='black', step_increase=0.1) +
  theme_minimal() +
  labs(x='Model', y='MCC') +
  theme(axis.text.x=element_blank()) +
  scale_colour_manual(values=model.colours, name='Model')
 

pdf('compare_HVG_across_models_celltypecol_update.pdf', width=12)
#pdf('compare_HVG_across_models_celltypecol.pdf', width=12)
average_model_metrics_df = average_model_metrics_list$combined   %>%  subset(gene.set == 'HVG')
ggplot(average_model_metrics_df  , aes(x=model, y=MCC, colour=celltype)) +
  geom_jitter(width=0.2, size=3) +
  geom_boxplot(outlier.shape = NA, color = 'black', fill = NA) +
  theme_minimal() +
  labs(x='Model', y='MCC') + 
  scale_colour_manual(values=celltypes.colours, name='Celltype')
dev.off() 


 


pdf('compare_across_models_and_geneset.pdf')
average_model_metrics_df = average_model_metrics_list$combined 
#%>% subset(gene.set %in% c('chrX','HVG', 'SLE'))

ggplot(average_model_metrics_df  , aes(x=gene.set, y=MCC, colour=model)) +
  geom_jitter(width=0.2, size=3) +
  geom_boxplot(outlier.shape = NA, color = 'black', fill = NA) +
  theme_minimal() +
  labs(x='Gene set', y='MCC') +
  #theme(axis.text.x=element_blank()) +
  scale_colour_manual(values=model.colours, name='Model')
dev.off() 









### Read in metrics file for ensemble models across gene sets and splits ###
ensemble_metrics_list <- list()
for(i in 1:10){
  metric.files <- unlist(lapply(methods, function(method){
  list.files(paste0('perez_update/split_', i, '/', method, '/ensemble'), pattern='metrics_', full.names=TRUE)}))
  #list.files(paste0('pseudobulk/split_', i, '/', method, '/ensemble'), pattern='metrics_', full.names=TRUE)}))

  metric.files <- metric.files[file.size(metric.files) > 0]
  metrics <- lapply(metric.files, read.csv)
  names(metrics) <- gsub('split_\\d+/|ensemble/metrics_|.csv', '', metric.files)  %>% 
  gsub('_combined', '', .) %>%
  gsub('/', '_', .)
  metrics_df <- bind_rows(metrics, .id='celltype') %>%
    mutate(
      gene.set=str_extract(celltype, "HVG.autosome|chrX|autosome|HVG|SLE"),
      method = str_extract(celltype, "boruta|enet|intersection|combined"),
      celltype=gsub('perez_update_', '', celltype),
      #celltype=gsub('pseudobulk_', '', celltype),
      celltype=gsub('combined_|.HVG.autosome', '', celltype),
      celltype = gsub("^(boruta_|enet_|intersection_|combined_)", "", celltype),
      celltype=gsub('.SLE|.chrX|.autosome|.HVG', '', celltype),
      celltype = gsub("_", " ", celltype))
  # Filter out autosome, HVG, and HVG.autosome gene sets
  #metrics_df <- subset(metrics_df, gene.set %in% c('chrX', 'SLE'))
  ensemble_metrics_list[[i]] <- metrics_df
}
names(ensemble_metrics_list) <- paste('split', 1:10, sep='_')

 

# Combine metrics from all splits
ensemble_metrics_df <- bind_rows(ensemble_metrics_list, .id='split')
# Format gene.set as factor
ensemble_metrics_df$gene.set <- factor(ensemble_metrics_df$gene.set)
ensemble_metrics_df$celltype <- factor(ensemble_metrics_df$celltype, levels = celltypes)
# Format split as factor
ensemble_metrics_df$split <- factor(ensemble_metrics_df$split, levels=paste('split', 1:10, sep='_'))
# Add ensemble as method
ensemble_metrics_df$model <- 'ensemble'


average_ensemble_model_metrics_list = list() 
for (methodi in methods){  
  average_ensemble_model_metrics_df <- ensemble_metrics_df %>% 
    subset(method == methodi)  %>%
    group_by(celltype, gene.set,model) %>%
    summarise(
      Accuracy = mean(Accuracy),
      Precision = mean(Precision),
      Recall= mean(Recall),
      F1=mean(F1),
      F1_lower=mean(F1_lower),
      F1_upper=mean(F1_upper),
      AUC=mean(AUC),
      AUC_lower=mean(AUC_lower),
      AUC_upper=mean(AUC_upper),
      AUPRC=mean(AUPRC),
      AUPRC_lower=mean(AUPRC_lower),
      AUPRC_upper=mean(AUPRC_upper),
      Kappa=mean(Kappa),
      MCC=mean(MCC),
      MCC_lower=mean(MCC_lower),
      MCC_upper=mean(MCC_upper),
      n_features=round(mean(n_features),0)) %>%
    data.frame()
  
  average_ensemble_model_metrics_df$gene.set <- factor(average_ensemble_model_metrics_df$gene.set)
  average_ensemble_model_metrics_df$celltype <- factor(average_ensemble_model_metrics_df$celltype)
  average_ensemble_model_metrics_list[[methodi]] = average_ensemble_model_metrics_df
}
  
 



average_ensemble_model_metrics_df = average_ensemble_model_metrics_list$combined
# %>%  subset(gene.set == 'SLE')
   ggplot(average_ensemble_model_metrics_df, aes(x=gene.set, y=MCC, colour=model)) +
    geom_jitter(width=0.2, size=3) +
    geom_boxplot(outlier.shape = NA, colour='black', fill=NA) +
    theme_minimal() +
    labs(x='Gene Set', y='MCC') + 
    scale_colour_manual(values=gene.set.colours, name='Gene Set') 
  
  
  
  
 
pdf("compare_avg_mcc_combined_hvg_fig2_panelA.pdf", width=12)   
   average_model_metrics_df = average_model_metrics_list$combined    
   average_ensemble_model_metrics_df = average_ensemble_model_metrics_list$combined  
   average_ensemble_model_metrics_df_combined = rbind(average_model_metrics_df,average_ensemble_model_metrics_df )      
  
   average_ensemble_model_metrics_df_combined_temp = average_ensemble_model_metrics_df_combined  %>%  subset(gene.set == 'HVG')
   ggplot(average_ensemble_model_metrics_df_combined_temp  , aes(x=model, y=MCC, colour=celltype)) +
     geom_jitter(width=0.2, size=3) +
     geom_boxplot(outlier.shape = NA, color = 'black', fill = NA) +
     theme_minimal() +
     labs(x='Model', y='MCC') + 
     scale_colour_manual(values=celltypes.colours, name='Cell type')
dev.off()   



average_model_metrics_df = average_model_metrics_list$combined    
average_ensemble_model_metrics_df = average_ensemble_model_metrics_list$combined  
average_ensemble_model_metrics_df_combined = rbind(average_model_metrics_df,average_ensemble_model_metrics_df )      

average_ensemble_model_metrics_df_combined_temp = average_ensemble_model_metrics_df_combined   %>%  subset(gene.set == 'chrX')
p1 <- ggplot(average_ensemble_model_metrics_df_combined_temp  , aes(x=model, y=MCC, colour=celltype)) +
  geom_jitter(width=0.2, size=3) + ggtitle("chrX") + 
  geom_boxplot(outlier.shape = NA, color = 'black', fill = NA) +
  theme_minimal() +   theme(legend.position = "none") + 
  labs(x='Model', y='MCC') + ylim(0,1)  + 
  scale_colour_manual(values=celltypes.colours, name='Cell type')
average_ensemble_model_metrics_df_combined_temp = average_ensemble_model_metrics_df_combined   %>%  subset(gene.set == 'SLE')
p2 <- ggplot(average_ensemble_model_metrics_df_combined_temp  , aes(x=model, y=MCC, colour=celltype)) +
  geom_jitter(width=0.2, size=3) +ggtitle("SLE") + 
  geom_boxplot(outlier.shape = NA, color = 'black', fill = NA) +
  theme_minimal() +   theme(legend.position = "none") + 
  labs(x='Model', y='MCC') + ylim(0,1)  + 
  scale_colour_manual(values=celltypes.colours, name='Cell type')



average_ensemble_model_metrics_df_combined_temp = average_ensemble_model_metrics_df_combined   %>%  subset(gene.set == 'HVG')
p3 <- ggplot(average_ensemble_model_metrics_df_combined_temp  , aes(x=model, y=MCC, colour=celltype)) +
  geom_jitter(width=0.2, size=3) +ggtitle("HVG") + 
  geom_boxplot(outlier.shape = NA, color = 'black', fill = NA) +
  theme_minimal() +   theme(legend.position = "none") + 
  labs(x='Model', y='MCC') + ylim(0,1)  + 
  scale_colour_manual(values=celltypes.colours, name='Cell type')


pdf("compare_avg_mcc_combined_hvg_fig2_panelA_v2.pdf", width=12)   
average_ensemble_model_metrics_df_combined_temp = average_ensemble_model_metrics_df_combined   %>%  subset(gene.set == 'HVG')
 ggplot(average_ensemble_model_metrics_df_combined_temp  , aes(x=model, y=MCC, colour=celltype)) +
  geom_jitter(width=0.2, size=3) +ggtitle("HVG") + 
  geom_boxplot(outlier.shape = NA, color = 'black', fill = NA) +
  theme_minimal() +   
  labs(x='Model', y='MCC') + ylim(0,1)  + 
  scale_colour_manual(values=celltypes.colours, name='Cell type')
dev.off() 

##########################################################################
# Figure 2 
# Panel C 
pdf("compare_avg_mcc_combined_X_SLE_fig2_panelC_updated.pdf")
temp = average_ensemble_model_metrics_list$combined %>%  subset(gene.set %in% c('chrX', 'SLE'))
ggplot(temp, aes(x=MCC, y=celltype)) +
  geom_errorbarh(aes(xmin = MCC_lower, xmax = MCC_upper)) +
  geom_vline(xintercept = 0.7, linetype = 'dotted', color='red') +
  geom_vline(xintercept = 0.5, linetype = 'dotted', color='white') +
  geom_point() +
  # geom_point(data=best_model_point, color='red') +
  labs(x='MCC', y='') +
  theme(axis.text.y=element_text(size=12)) +
  facet_wrap(~gene.set, ncol=5, nrow=1)
 dev.off() 
 
 
 ensemble_metrics_df_sub = ensemble_metrics_df %>% subset(gene.set %in% c('chrX', 'SLE')) %>% subset(method == "combined") 
 compare_celltype <- lapply(split(ensemble_metrics_df_sub, ensemble_metrics_df_sub$celltype), function(x){
   wilcox.test(MCC ~ gene.set, data=x)$p.value
 })
 
 compare_celltype <- t(bind_rows(compare_celltype, .id='celltype'))
 compare_celltype <- data.frame(p.value=compare_celltype[,1], FDR=p.adjust(compare_celltype[,1], method='fdr'))
 subset(compare_celltype, FDR > 0.05)
 top_celltypes <- rownames(subset(compare_celltype, FDR > 0.05))
  
 
 top_celltypes_top_features <- selected_features$top_features[names(selected_features$top_features) %in% paste0(gsub(' ', '_', top_celltypes), '.chrX')]
 top_celltypes_top_features_mtx <- fromList(top_celltypes_top_features)
 rownames(top_celltypes_top_features_mtx) <- unique(unlist(top_celltypes_top_features))
 colnames(top_celltypes_top_features_mtx) <- gsub('.chrX', '', names(top_celltypes_top_features_mtx)) %>% gsub('_', ' ', .)
 
 
 celltypes_average <- subset(average_ensemble_model_metrics_list$combined , gene.set == 'SLE')[,c('celltype', 'MCC', 'MCC_lower', 'MCC_upper')]
 celltypes_average2 <- subset(average_ensemble_model_metrics_list$combined , gene.set == 'chrX')[,c('celltype', 'MCC', 'MCC_lower', 'MCC_upper')]
 
 top_celltypes_average <- subset(average_ensemble_model_metrics_list$combined , celltype %in% top_celltypes & gene.set == 'chrX')[,c('celltype', 'MCC', 'MCC_lower', 'MCC_upper')]
 # round 'MCC', 'MCC_lower', 'MCC_upper' to 2 decimal
 top_celltypes_average$MCC <- round(top_celltypes_average$MCC, 2)
 top_celltypes_average$MCC_lower <- round(top_celltypes_average$MCC_lower, 2)
 top_celltypes_average$MCC_upper <- round(top_celltypes_average$MCC_upper, 2)
 
 write.csv(top_celltypes_average, 'top_celltypes_average_updated.csv')
 
 subset(average_ensemble_model_metrics_list$combined , gene.set == 'chrX' & MCC_upper > 0.7 & celltype %in% top_celltypes)[,c('celltype', 'MCC', 'MCC_lower', 'MCC_upper')]
 
 
 lapply(split(ensemble_metrics_df_sub, ensemble_metrics_df_sub$celltype), function(x){
   pairwise.wilcox.test(x$MCC, x$gene.set, p.adjust.method='fdr')$p.value
 })
 
 # Calculate p-values for each cell type comparison (all celltypes)
 pval_list <- lapply(split(ensemble_metrics_df_sub, ensemble_metrics_df_sub$celltype), function(x){
   wilcox.test(MCC ~ gene.set, data=x)$p.value
 })
 
 
 all_celltype_pvals <- data.frame(
   celltype = names(pval_list),
   p_value = unlist(pval_list)
 ) %>%
   mutate(p_signif = case_when(
     p_value < 0.001 ~ "***",
     p_value < 0.01 ~ "**", 
     p_value < 0.05 ~ "*",
     TRUE ~ "ns"
   ))
 
 
 comparisons <- list(c('chrX', 'SLE'))

 pdf("compare_avg_mcc_combined_X_SLE_fig2_panelC_v2.pdf")
 
  ggplot(ensemble_metrics_df_sub, aes(x=MCC, y=celltype, fill=gene.set)) +
   geom_boxplot(position=position_dodge(width=0.8), alpha=0.7) +
   geom_point(position=position_jitterdodge(dodge.width=0.8, jitter.width=0.2), 
              size=1) +
   geom_text(data=all_celltype_pvals, aes(x=Inf, y=celltype, label=p_signif), 
             hjust=1.2, vjust=0.5, size=3, inherit.aes=FALSE) +
  geom_vline(xintercept = 0.7, linetype = 'dotted', color='red') + 
   theme_minimal() +
   labs(x='MCC', y='') +
   theme(
     axis.text.y=element_text(size=14),
     axis.text.x=element_text(size=10 ),
     legend.position="top") +
   scale_fill_manual(values=gene.set.colours, name='Gene Set')

 dev.off()   
 
 
 library(forcats)
 o2 = order(celltypes_average[,4])
 o3 = order(celltypes_average2[,2]) 
 o4 = order(all_celltype_pvals$p_value) 
 
 temp = ensemble_metrics_df_sub 
 temp$celltype=  fct_relevel(ensemble_metrics_df_sub$celltype , 
                             celltypes[rev(o4)] ) 
    
 #
 pdf("compare_avg_mcc_combined_X_SLE_fig2_panelC_v3.pdf")
 pdf("compare_avg_mcc_combined_X_SLE_fig2_panelC_updated_v4.pdf")
   ggplot(temp, aes(y=MCC, x=celltype, fill=gene.set)) +
   geom_boxplot(position=position_dodge(width=0.8), alpha=0.7) +
   geom_point(position=position_jitterdodge(dodge.width=0.8, jitter.width=0.2), 
              size=1) +
   geom_text(data=all_celltype_pvals, aes(y=Inf, x=celltype, label=p_signif), 
             hjust=1.2, vjust=0.5, size=3, inherit.aes=FALSE) +
   geom_hline(yintercept = 0.7, linetype = 'dotted', color='red') + 
   theme_minimal() +
   labs(y='MCC', x='') +
   theme(
     axis.text.y=element_text(size=10),
     axis.text.x=element_text(size=14, angle=90),
     legend.position="top") +
   scale_fill_manual(values=gene.set.colours, name='Gene Set')
 
 dev.off()   
 
 
 
 
##########################################################################

# Read in features for all ML models across gene sets and splits
models_features_list <- list()
for(i in 1:10){
  features.files <- unlist(lapply(methods, function(method){
    list.files(paste0('perez_update/split_', i, '/features/'), pattern=method, full.names=TRUE) }))
  #  list.files(paste0('pseudobulk/split_', i, '/features/'), pattern=method, full.names=TRUE) }))

  # Check if metric files are empty
  features.files <- features.files[file.size(features.files) > 0]
  features <- lapply(features.files, read.csv)
  names(features) <- gsub('split_\\d+/|metrics/|.csv', '', features.files) %>% 
    gsub('_metrics', '', .) %>%
    gsub('/', '_', .)
  features_df <- bind_rows(features, .id = 'celltype') %>%
    mutate(
      gene.set = str_extract(celltype, "HVG\\.autosome|chrX|autosome|HVG|SLE"),
      method = str_extract(celltype, "boruta|enet|intersection|combined"),
      celltype = gsub("^perez_update_features_", "", celltype),
      #celltype = gsub("^pseudobulk_features_", "", celltype),
      
      celltype = gsub("features.", "", celltype),
      celltype = gsub("(boruta_|enet_|intersection_|combined_)", "", celltype),
      celltype = gsub("\\.HVG\\.autosome|\\.SLE|\\.chrX|\\.autosome|\\.HVG", "", celltype),
      celltype = gsub("_", " ", celltype)
    )
 
  models_features_list[[i]] <- features_df
}

names(models_features_list) <- paste('split', 1:10, sep='_')


# Combine metrics from all splits
features_list_df <- bind_rows(models_features_list, .id='split')
# Format gene.set as factor
features_list_df$gene.set <- factor(features_list_df$gene.set)
features_list_df$celltype <- factor(features_list_df$celltype, levels = celltypes)
# Format split as factor
features_list_df$split <- factor(features_list_df$split, levels=paste('split', 1:10, sep='_'))
 



features_list_df_temp = features_list_df %>%  subset(method == 'combined')
features_temp <-  lapply(celltypes, function(celltypei) 
  features_list_df %>%  subset(method == 'combined') %>%  subset(celltype == celltypei)) 
names(features_temp) = celltypes 
 
features_temp2 <-  lapply(celltypes, function(celltypei) 
    tapply(features_temp[[celltypei]]$Feature, features_temp[[celltypei]]$gene.set, plyr::count))
names(features_temp2) = celltypes 


all_features = unique(unlist(lapply(celltypes, function(i)  features_temp2[[i]]$chrX[,1])))
top_features = unique(unlist(lapply(celltypes, function(i)   features_temp2[[i]]$chrX[features_temp2[[i]]$chrX[,2]==10,1]   ))) 

temp = do.call(rbind, lapply(celltypes, function(i)  cbind(i, features_temp2[[i]]$chrX))  ) 

all_celltypes_all_features_mtx = spread( temp, key=1, value = 3, fill = 0 )
rownames(all_celltypes_all_features_mtx) = all_celltypes_all_features_mtx[,1] 
temp_mtx = apply(all_celltypes_all_features_mtx[,-1], 2, as.numeric)
rownames(temp_mtx ) = rownames(all_celltypes_all_features_mtx)
all_celltypes_all_features_mtx = temp_mtx 

heatmap.2(all_celltypes_all_features_mtx , density.info = "none", trace = "none", col=magma(10))


all_celltypes_top_features_mtx = 1*(all_celltypes_all_features_mtx == 10 )


x_escape = rownames(escape)[escape$Xie_Tukiainen2017==1]
x_escape2 = rownames(escape)[escape$Xive_Tukiainen2017==1]

m = match(celltypes, colnames(all_celltypes_top_features_mtx))
all_celltypes_top_features_mtx = all_celltypes_top_features_mtx[,m]  
all_celltypes_all_features_mtx = all_celltypes_all_features_mtx[,m]  

xi_cols = rep("grey", length(rownames(all_celltypes_top_features_mtx) ))
xi_cols[!is.na(match(rownames(all_celltypes_top_features_mtx), x_escape)) ] = "purple"
xi_cols[!is.na(match(rownames(all_celltypes_top_features_mtx), x_escape2)) ] = "darkorchid4"
 
f1 = rowSums(all_celltypes_top_features_mtx) > 0 

f2 =  (pval_list > 0.05)

pdf("all_celltypes_top_features_heatmap_fig3_panelA.pdf")  
heatmap.2(all_celltypes_top_features_mtx[f1,] , 
          density.info = "none", trace = "none", 
          col=c(0,1),
          RowSideColors = xi_cols[f1],
          ColSideColors=celltypes.colours[])
dev.off() 

cols_grey = colorRampPalette(c("white", "black"))( 10 ) 

pdf("all_celltypes_all_features_heatmap_fig3_panelA.pdf")  
heatmap.2(all_celltypes_all_features_mtx[f1,] , 
          density.info = "none", trace = "none", 
          col= cols_grey,
          RowSideColors = xi_cols[f1],
          ColSideColors=celltypes.colours[])
dev.off() 




heatmap.2(all_celltypes_top_features_mtx[f1,f2] , 
          density.info = "none", trace = "none", 
          col=c(0,1),
          RowSideColors = xi_cols[f1],
          ColSideColors=celltypes.colours[f2])




heatmap.2(all_celltypes_all_features_mtx[f1,f2] , 
          density.info = "none", trace = "none", 
          col=inferno(10),
          RowSideColors = xi_cols[f1],
          ColSideColors=celltypes.colours[f2])



heatmap.2(all_celltypes_all_features_mtx , 
          density.info = "none", trace = "none", 
          col=inferno(10),
          RowSideColors = xi_cols,
          ColSideColors=celltypes.colours)



heatmap.2(all_celltypes_all_features_mtx[,f2] , 
          density.info = "none", trace = "none", 
          col=inferno(10),
          RowSideColors = xi_cols[],
          ColSideColors=celltypes.colours[f2])






########################################### 
load("edgeR.all_celltypes.RData")
names(edgeR) = gsub("_", " ", names(edgeR)) 
thresh = 0.5 
degsn = sapply(celltypes, function(i)  sum(abs(edgeR[[i]]$logFC) > thresh & edgeR[[i]]$FDR <0.05) ) 
updegsn = sapply(celltypes, function(i)  sum( (edgeR[[i]]$logFC) > thresh & edgeR[[i]]$FDR <0.05) ) 
downdegsn = sapply(celltypes, function(i)  sum( (edgeR[[i]]$logFC) < -thresh & edgeR[[i]]$FDR <0.05) ) 

degs = lapply(celltypes, function(i)   edgeR[[i]][(abs(edgeR[[i]]$logFC) > thresh & edgeR[[i]]$FDR <0.05),1] ) 
updegs = lapply(celltypes, function(i)   edgeR[[i]][( (edgeR[[i]]$logFC) > thresh & edgeR[[i]]$FDR <0.05),1] ) 
downdegs = lapply(celltypes, function(i)   edgeR[[i]][( (edgeR[[i]]$logFC) < -thresh & edgeR[[i]]$FDR <0.05),1] ) 

names(degs) = celltypes 
names(updegs) = celltypes 
names(downdegs) = celltypes 



deg_table = cbind(degsn,updegsn, downdegsn )

o1 = order(deg_table[,1])

pdf("DEGS_total_fig2_panelB.pdf")
bp = barplot( t(deg_table[rev(o1),2:3]) , col=cividis(2), border = 0, hor=T) 
text(500, bp, rownames(deg_table)[rev(o1)], pos=2 )
dev.off() 



features_overlap <-  lapply(celltypes, function(celltypei) 
                      lapply(gene.set, function(gi)  
                        intersect(features_temp2[[celltypei]][[gi]][,1],  degs[[celltypei]])))


features_overlapn <-  lapply(celltypes, function(celltypei) 
  lapply(gene.set, function(gi)  
    length(intersect(features_temp2[[celltypei]][[gi]][,1],  degs[[celltypei]]))))


names(features_overlapn ) = celltypes
names(features_overlap ) = celltypes

 for( celltypei in celltypes){ 
   
   names(features_overlap[[celltypei]]) = gene.set
   names(features_overlapn[[celltypei]]) = gene.set
   } 
 

deg_feature_overlaps = t(sapply(celltypes, function(celltypei) unlist(features_overlapn[[celltypei]]) )) 


chrX_overlap = lapply(celltypes, function(celltypei) intersect(degs[[celltypei]], chrX))
escape_overlap = lapply(celltypes, function(celltypei) intersect(degs[[celltypei]], rownames(escape)[escape$Xie_Tukiainen2017==1] ))

names(chrX_overlap) = celltypes 
names(escape_overlap) = celltypes 



chrX_overlapn = sapply(celltypes, function(celltypei) length(intersect(degs[[celltypei]], chrX)))
escape_overlapn = sapply(celltypes, function(celltypei) length(intersect(degs[[celltypei]], rownames(escape)[escape$Xie_Tukiainen2017==1]) ))

cbind( deg_table, chrX_overlapn, escape_overlapn)



celltypei= "Classical monocyte"
oo = order(edgeR[[celltypei]]$logFC, -log10(edgeR[[celltypei]]$FDR))
x = 1:length(oo)
y2 = (edgeR[[celltypei]][oo,'logFC'] ) 
plot(1:length(oo), y2)
m = match(features_overlap[[celltypei]]$SLE , edgeR[[celltypei]][oo,1])
sle1 = m[!is.na(m)]
m = match(features_overlap[[celltypei]]$HVG , edgeR[[celltypei]][oo,1])
hvg1 = m[!is.na(m)]
m = match(features_overlap[[celltypei]]$chrX , edgeR[[celltypei]][oo,1])
chrX1 = m[!is.na(m)]
abline(v = sle1, col=gene.set.colours[['SLE']])
abline(v = hvg1, col=gene.set.colours[['HVG']])
abline(v = chrX1, col=gene.set.colours[['chrX']])



celltypei= "B cell"
#celltypei= "Classical monocyte"
#celltypei= "NK cell Bright"
celltypei= "CD8 positive, alpha beta T cell"

pdf("degs_overlaps_update.pdf")
#pdf("degs_overlaps.pdf")
for (celltypei in celltypes){ 
  oo = order(edgeR[[celltypei]]$logFC, -log10(edgeR[[celltypei]]$FDR))
  x = 1:length(oo)
  y2 = (edgeR[[celltypei]][oo,'logFC'] ) 
  plot(x, y2, col="grey", pch=19, cex=0.5, main=celltypei, xlab="Gene rank", ylab="log2FC")
  abline(h=c(0.5, -0.5), lty=2, col="grey" )
  m = match(features_overlap[[celltypei]]$SLE , edgeR[[celltypei]][oo,1])
  sle1 = m[!is.na(m)]
  m = match(features_overlap[[celltypei]]$HVG , edgeR[[celltypei]][oo,1])
  hvg1 = m[!is.na(m)]
  m = match(features_overlap[[celltypei]]$chrX , edgeR[[celltypei]][oo,1])
  chrX1 = m[!is.na(m)]
  
  if(length(chrX1) > 0 ) {
    abline(v = chrX1, col=gene.set.colours[['chrX']])
    points(x[chrX1], y2[chrX1], col=gene.set.colours[['chrX']], pch=19)
    text(x[chrX1], y2[chrX1], edgeR[[celltypei]][oo,][chrX1,1], pos=2, font=3, xpd=NA)
  }
  if(length(hvg1) > 0 ) { 
    abline(v = hvg1, col=gene.set.colours[['HVG']])
    points(x[hvg1], y2[hvg1], col=gene.set.colours[['HVG']], pch=19)
    text(x[hvg1], y2[hvg1], edgeR[[celltypei]][oo,][hvg1,1], pos=2, font=3, xpd=NA)
  }
  if(length(sle1) > 0 ) { 
    abline(v = sle1, col=gene.set.colours[['SLE']])
      points(x[sle1], y2[sle1], col=gene.set.colours[['SLE']], pch=19)
    text(x[sle1], y2[sle1], edgeR[[celltypei]][oo,][sle1,1], pos=2, font=3, xpd=NA)
  }
} 
dev.off()  



 deg_feature_overlaps_expanded = list() 
 
for (celltypei in celltypes){ 
  deg_feature_overlaps_expanded[[celltypei]] = list() 
  oo = order(edgeR[[celltypei]]$logFC, -log10(edgeR[[celltypei]]$FDR))
  m = match(features_overlap[[celltypei]]$SLE , edgeR[[celltypei]][oo,1])
  sle1 = m[!is.na(m)]
  m = match(features_overlap[[celltypei]]$HVG , edgeR[[celltypei]][oo,1])
  hvg1 = m[!is.na(m)]
  m = match(features_overlap[[celltypei]]$chrX , edgeR[[celltypei]][oo,1])
  chrX1 = m[!is.na(m)]
  
  if(length(chrX1) > 0 ) {
    deg_feature_overlaps_expanded[[celltypei]][['chrX']] = cbind(edgeR[[celltypei]][oo,][chrX1,],oo[chrX1])
  }
  if(length(hvg1) > 0 ) { 
    deg_feature_overlaps_expanded[[celltypei]][['HVG']] =  cbind(edgeR[[celltypei]][oo,][hvg1,],oo[hvg1])
  }
  if(length(sle1) > 0 ) { 
    deg_feature_overlaps_expanded[[celltypei]][['SLE']] =  cbind(edgeR[[celltypei]][oo,][sle1,], oo[sle1])
  }
}  


 for(celltypei in celltypes) { temp = features_list_df %>% subset(method == "combined") %>% subset( gene.set == "HVG") %>% subset( celltype == celltypei); ft = grep("[a-z]", temp$Feature , invert=T); temp2 = plyr::count(temp[ft,1]) ; print(c(celltypei,mean(temp2[,2])) )  } 
 for(celltypei in celltypes) { temp = features_list_df %>% subset(method == "combined") %>% subset( gene.set == "SLE") %>% subset( celltype == celltypei); ft = grep("[a-z]", temp$Feature , invert=T); temp2 = plyr::count(temp[ft,1]) ; print(c(celltypei,mean(temp2[,2])) )  } 
 for(celltypei in celltypes) { temp = features_list_df %>% subset(method == "combined") %>% subset( gene.set == "chrX") %>% subset( celltype == celltypei); ft = grep("[a-z]", temp$Feature , invert=T); temp2 = plyr::count(temp[ft,1]) ; print(c(celltypei,mean(temp2[,2])) )  } 

 
 
sapply(celltypes, function(celltypei) length(edgeR[[celltypei]][,1]))

 

################################### 
load("phens.Rdata")

phens.df = as.data.frame(phens)
colnames(phens.df) = c("donorID", "age", "sex", "ancestry", "status")
phens.df$age = as.numeric(phens.df$age)
par(mfrow=c(2,2))
i = 1 
for( status in unique(phens.df$status)){ 
  for( ancestry in unique(phens.df$ancestry)){ 
    xtemp = phens.df$age[phens.df$ancestry == ancestry & phens.df$status == status] 
    hist(xtemp, breaks=30, xlim=c(20,80), ylim=c(0,15), col=viridis(5)[i], border=0,
         main=paste(ancestry, status), xlab="Age")
    abline(v = mean(xtemp), lty=2, col=viridis(5)[i] )
    i = i + 1 
  }
}
 


freq_d  = plyr:: count( metadata$cell_type_detailed2[metadata$disease=="disease" & !is.na(metadata$cell_type_detailed2)])
freq_c  = plyr:: count( metadata$cell_type_detailed2[metadata$disease=="control" & !is.na(metadata$cell_type_detailed2)] )
freq    = plyr:: count( metadata$cell_type_detailed2[!is.na(metadata$cell_type_detailed2)] )
freq_a  = plyr:: count( metadata$cell_type_detailed2[metadata$ancestry=="Asian"    & !is.na(metadata$cell_type_detailed2)] )
freq_e  = plyr:: count( metadata$cell_type_detailed2[metadata$ancestry=="European" & !is.na(metadata$cell_type_detailed2)] )
freq_ac = plyr:: count( metadata$cell_type_detailed2[metadata$ancestry=="Asian"    & metadata$disease=="control" & !is.na(metadata$cell_type_detailed2)] )
freq_ad = plyr:: count( metadata$cell_type_detailed2[metadata$ancestry=="Asian"    & metadata$disease=="disease" & !is.na(metadata$cell_type_detailed2)] )
freq_ec = plyr:: count( metadata$cell_type_detailed2[metadata$ancestry=="European" & metadata$disease=="control" & !is.na(metadata$cell_type_detailed2)] )
freq_ed = plyr:: count( metadata$cell_type_detailed2[metadata$ancestry=="European" & metadata$disease=="disease" & !is.na(metadata$cell_type_detailed2)] )

freq_stat_cell = cbind(freq, freq_d[,2], freq_c[,2],  freq_a[,2], freq_e[,2],freq_ad[,2], freq_ed[,2] , freq_ac[,2], freq_ec[,2]  )
colnames(freq_stat_cell) = c("celltype", "total","disease", "control",  "Asian", "European",
                             "disease - Asian","disease - European", 
                             "control - Asian","control - European")

save(freq_stat_cell, file="freq_stat_cell.Rdata")




m = match( names(celltypes.colours), freq_stat_cell[,1])
freq_stat_cell = freq_stat_cell[m,]
barplot( log10( t(t(freq_stat_cell[,-1] )) ) , col=celltypes.colours , beside=T, hor=T ) 
 



filt = !is.na(metadata$cell_type_detailed2)  
freq1 = plyr::count( cbind( metadata$donor_id ,as.character(metadata$cell_type_detailed2))[filt,] )

freq_p0 = spread(freq1, key = "x.2", value = "freq")
tenmp = t(freq_p0[ ,-1])
tenmp[is.na(tenmp)] = 0
tenmp = as.matrix(tenmp)
tempfrac0 = (sapply(1:dim(tenmp)[2], function(i) tenmp[,i]/sum(tenmp[,i], na.rm=T) ) )
colnames(tempfrac0) = freq_p0[,1]
m = match(colnames(tempfrac0), phens[,1] )
phens = phens[m,]

props0 =  propeller(clusters=metadata$cell_type_detailed2[filt], sample=metadata$donor_id[filt], group=paste(metadata$disease[filt],metadata$ancestry[filt]))

m = match( rownames(tempfrac0),rownames(props0) )
props0 = props0[m,]

tempfrac0l_ec = tempfrac0[,phens[,4]=="European" & phens[,5]=="control"]
tempfrac0l_ac = tempfrac0[,phens[,4]=="Asian" & phens[,5]=="control" ] 
tempfrac0l_ed = tempfrac0[,phens[,4]=="European" & phens[,5]=="disease"]
tempfrac0l_ad = tempfrac0[,phens[,4]=="Asian" & phens[,5]=="disease" ] 



bw=0.05
tempfrac0l_ec.d = lapply(1:dim(tempfrac0)[1], function(i) density(tempfrac0[i, phens[,4]=="European" & phens[,5]=="control"], bw=bw, from=0, to=1) )
tempfrac0l_ac.d = lapply(1:dim(tempfrac0)[1], function(i) density(tempfrac0[i, phens[,4]=="Asian" & phens[,5]=="control"], bw=bw, from=0, to=1) )
tempfrac0l_ed.d = lapply(1:dim(tempfrac0)[1], function(i) density(tempfrac0[i, phens[,4]=="European" & phens[,5]=="disease"], bw=bw, from=0, to=1) )
tempfrac0l_ad.d = lapply(1:dim(tempfrac0)[1], function(i) density(tempfrac0[i, phens[,4]=="Asian" & phens[,5]=="disease"], bw=bw, from=0, to=1) )
tempfrac0l.d    = lapply(1:dim(tempfrac0)[1], function(i) density(tempfrac0[i,], bw=bw, from=0, to=1) )

tempfrac0l.d2    = lapply(1:length(tempfrac0l.d),  function(i) cbind( c(0,tempfrac0l.d[[i]]$x,1),  c(0,tempfrac0l.d[[i]]$y,0) ) )
tempfrac0l_ec.d2 = lapply(1:length(tempfrac0l_ec.d),  function(i) cbind( c(0,tempfrac0l_ec.d[[i]]$x,1),  c(0,tempfrac0l_ec.d[[i]]$y,0) ) )
tempfrac0l_ac.d2 = lapply(1:length(tempfrac0l_ac.d),  function(i) cbind( c(0,tempfrac0l_ac.d[[i]]$x,1),  c(0,tempfrac0l_ac.d[[i]]$y,0) ) )
tempfrac0l_ed.d2 = lapply(1:length(tempfrac0l_ed.d),  function(i) cbind( c(0,tempfrac0l_ed.d[[i]]$x,1),  c(0,tempfrac0l_ed.d[[i]]$y,0) ) )
tempfrac0l_ad.d2 = lapply(1:length(tempfrac0l_ad.d),  function(i) cbind( c(0,tempfrac0l_ad.d[[i]]$x,1),  c(0,tempfrac0l_ad.d[[i]]$y,0) ) )


save(tempfrac0,freq_p0, props0, phens, tempfrac0l.d, tempfrac0l.d2,
     tempfrac0l_ec, tempfrac0l_ac,tempfrac0l_ed, tempfrac0l_ad,
     tempfrac0l_ec.d, tempfrac0l_ac.d,tempfrac0l_ed.d, tempfrac0l_ad.d,
     tempfrac0l_ec.d2, tempfrac0l_ac.d2,tempfrac0l_ed.d2, tempfrac0l_ad.d2,
     file="props_tests_cell.Rdata")


 

nn = 75
n =  dim(tempfrac0)[1]  
ni = nn/n
nni = 5
oi = 1:n 

m = match( names(celltypes.colours), rownames(tempfrac0) )
f.c = !is.na(m)
f.t = m[f.c] 

celltypes.colours2 = celltypes.colours[f.c]


tempfrac0   =   tempfrac0[f.t,]
tempfrac0l_ac =   tempfrac0l_ac[f.t,]
tempfrac0l_ec =   tempfrac0l_ec[f.t,]
tempfrac0l_ad =   tempfrac0l_ad[f.t,]
tempfrac0l_ed =   tempfrac0l_ed[f.t,]



tempfrac0l_ac.d   =   tempfrac0l_ac.d[f.t]
tempfrac0l_ac.d2  =   tempfrac0l_ac.d2[f.t]
tempfrac0l_ec.d   =   tempfrac0l_ec.d[f.t]
tempfrac0l_ec.d2  =   tempfrac0l_ec.d2[f.t]

tempfrac0l_ad.d   =   tempfrac0l_ad.d[f.t]
tempfrac0l_ad.d2  =   tempfrac0l_ad.d2[f.t]
tempfrac0l_ed.d   =   tempfrac0l_ed.d[f.t]
tempfrac0l_ed.d2  =   tempfrac0l_ed.d2[f.t]
 
props0 = props0[f.t,]

xxx = 0*(props0$FDR)
xxx[(props0$FDR  < 0.05)] = "*"
xxx[(props0$FDR  < 0.01)] = "**"
xxx[(props0$FDR  < 0.001)] = "***"
xxx[xxx==0] = ""




plot(-10,-10, xlim=c(0,1), ylim=c(0,nn+ni), xlab="Fraction", ylab="Density" ,axes=F)
axis(1)
tt = lapply(oi, function(i) abline( h = ((i-1)*ni)  ,  col="grey", lty=2 ) ) 
tt = lapply(oi, function(i) polygon(tempfrac0l_ac.d2[[i]][,1],((i-1)*ni)+(tempfrac0l_ac.d2[[i]][,2]*nni/n) ,  border=NA, col= (celltypes.colours2[i])))
tt = lapply(oi, function(i) polygon(tempfrac0l_ec.d2[[i]][,1],((i-1)*ni)+(tempfrac0l_ec.d2[[i]][,2]*nni/n) ,  border=NA, col= (celltypes.colours2[i])))
tt = lapply(oi, function(i) polygon(tempfrac0l_ad.d2[[i]][,1],((i-1)*ni)+(tempfrac0l_ad.d2[[i]][,2]*nni/n) ,  border=NA, col= (celltypes.colours2[i])))
tt = lapply(oi, function(i) polygon(tempfrac0l_ed.d2[[i]][,1],((i-1)*ni)+(tempfrac0l_ed.d2[[i]][,2]*nni/n) ,  border=NA, col= (celltypes.colours2[i])))
tt = lapply(oi, function(i) segments( mean(tempfrac0l_ac[i,]),  ((i-1)*ni) , 
                                      mean(tempfrac0l_ac[i,]), ((i-1)*ni) + ni/2 , lwd=2, col= viridis(5)[1]))
tt = lapply(oi, function(i) segments( mean(tempfrac0l_ec[i,]),  ((i-1)*ni) , 
                                      mean(tempfrac0l_ec[i,]), ((i-1)*ni) + ni/2 , lwd=2, col= viridis(5)[2]))
tt = lapply(oi, function(i) segments( mean(tempfrac0l_ad[i,]),  ((i-1)*ni) , 
                                      mean(tempfrac0l_ad[i,]), ((i-1)*ni) + ni/2 , lwd=2, col= viridis(5)[3]))
tt = lapply(oi, function(i) segments( mean(tempfrac0l_ed[i,]),  ((i-1)*ni) , 
                                      mean(tempfrac0l_ed[i,]), ((i-1)*ni) + ni/2 , lwd=2, col= viridis(5)[4]))
text(0.8, ((oi -1)*ni), names(celltypes.colours2))
text(1, ((oi -1)*ni), xxx)



pdf("props_tests_cell.pdf")
xmax =  (max(tempfrac0))
plot(-10,-10, xlim=c(0,1), ylim=c(0,nn+ni), xlab="Fraction", ylab="Density" ,axes=F)
axis(1)
tt = lapply(oi, function(i) abline( h = ((i-1)*ni)  ,  col="grey", lty=2 ) ) 
tt = lapply(oi, function(i) polygon(tempfrac0l_ac.d2[[i]][,1],((i-1)*ni)+(tempfrac0l_ac.d2[[i]][,2]*nni/n) ,  border=NA, col= (celltypes.colours2[i])))
tt = lapply(oi, function(i) polygon(tempfrac0l_ec.d2[[i]][,1],((i-1)*ni)+(tempfrac0l_ec.d2[[i]][,2]*nni/n) ,  border=NA, col= (celltypes.colours2[i])))
tt = lapply(oi, function(i) polygon(tempfrac0l_ad.d2[[i]][,1],((i-1)*ni)+(tempfrac0l_ad.d2[[i]][,2]*nni/n) ,  border=NA, col= (celltypes.colours2[i])))
tt = lapply(oi, function(i) polygon(tempfrac0l_ed.d2[[i]][,1],((i-1)*ni)+(tempfrac0l_ed.d2[[i]][,2]*nni/n) ,  border=NA, col= (celltypes.colours2[i])))
tt = lapply(oi, function(i) segments( mean(tempfrac0l_ac[i,]),  ((i-1)*ni) , 
                                      mean(tempfrac0l_ac[i,]), ((i-1)*ni) + ni/2 , lwd=2, col= viridis(5)[1]))
tt = lapply(oi, function(i) segments( mean(tempfrac0l_ec[i,]),  ((i-1)*ni) , 
                                      mean(tempfrac0l_ec[i,]), ((i-1)*ni) + ni/2 , lwd=2, col= viridis(5)[2]))
tt = lapply(oi, function(i) segments( mean(tempfrac0l_ad[i,]),  ((i-1)*ni) , 
                                      mean(tempfrac0l_ad[i,]), ((i-1)*ni) + ni/2 , lwd=2, col= viridis(5)[3]))
tt = lapply(oi, function(i) segments( mean(tempfrac0l_ed[i,]),  ((i-1)*ni) , 
                                      mean(tempfrac0l_ed[i,]), ((i-1)*ni) + ni/2 , lwd=2, col= viridis(5)[4]))
text(0.95, ((oi -1)*ni), names(celltypes.colours2), pos=2)
text(1, ((oi -1)*ni), xxx)
dev.off() 

phens = as.data.frame(phens)
f_ec = phens[,4]=="European" & phens[,5]=="control"
f_ac = phens[,4]=="Asian" & phens[,5]=="control" 
f_ed = phens[,4]=="European" & phens[,5]=="disease"
f_ad = phens[,4]=="Asian" & phens[,5]=="disease" 

phens[,2] =  as.numeric(phens[,2] )
nC = c(5,6)
n = dim(tempfrac0)[1] 
colnames(phens) = c("sample_uuid", "age", "sex", "ancestry", "disease") 


cors = sapply(1:n, function(i) cor.test( phens[,2], tempfrac0[i,]  , m= "s")  ) 
cors_e = unlist(cors[4,])
pvals = unlist(cors[3,] ) 
cors_all = cbind(rownames(tempfrac0) , pvals, p.adjust(pvals), cors_e)
collabels = cbind(rownames(tempfrac0) ,  celltypes.colours) 


cors_ac = sapply(1:n, function(i) cor.test( phens$age[f_ac], tempfrac0[i,f_ac]  , m= "s")  ) 
cors_ac_e = unlist(cors_ac[4,])
pvals_ac = unlist(cors_ac[3,] ) 

cors_ad = sapply(1:n, function(i) cor.test( phens$age[f_ad], tempfrac0[i,f_ad]  , m= "s")  ) 
cors_ad_e = unlist(cors_ad[4,])
pvals_ad = unlist(cors_ad[3,] ) 


cors_ec = sapply(1:n, function(i) cor.test( phens$age[f_ec], tempfrac0[i,f_ec]  , m= "s")  ) 
cors_ec_e = unlist(cors_ec[4,])
pvals_ec = unlist(cors_ec[3,] ) 

cors_ed = sapply(1:n, function(i) cor.test( phens$age[f_ed], tempfrac0[i,f_ed]  , m= "s")  ) 
cors_ed_e = unlist(cors_ed[4,])
pvals_ed = unlist(cors_ed[3,] ) 



cors_all = cbind( cors_all, 
                  pvals_ac, p.adjust(pvals_ac), cors_ac_e,
                  pvals_ec, p.adjust(pvals_ec), cors_ec_e,  
                  pvals_ad, p.adjust(pvals_ad), cors_ad_e,  
                  pvals_ed, p.adjust(pvals_ed), cors_ed_e)  
                  


temp = apply(cors_all[,2:16],2, as.numeric ) 
rownames(temp) = rownames(tempfrac0) 
cors_all = temp 
save(cors_all, file="cors_counts_age.Rdata") 

a =  (c(cors_all[,3], cors_all[,6], cors_all[,9], cors_all[,12], cors_all[,15]) )
b =  (c(rep("All",n), rep("Asian - control",n), rep("European - control",n), rep("Asian - disease",n), rep("European - disease",n) )) 
d = rep(celltypes.colours2,5) 


aa =  -log10(c(cors_all[,2], cors_all[,5], cors_all[,8], cors_all[,11], cors_all[,14])  )  
bb = 2 *(aa >= -log10(0.05) ) + 1 
beeswarm(a~b, pch=19, pwcol=d , ylab="Correlations (rho)", xlab="", bty="n",pwcex=bb) 

 







# Read in metrics for all ML models across gene sets and splits for validation set
models_metrics_list <- list()
metric.files  <- list.files(paste0('nehar_belaid_update/'), pattern='metrics_', full.names=TRUE) 
metric.files <- metric.files[file.size(metric.files) > 0]
metrics <- lapply(metric.files, read.csv)
names(metrics) <- gsub('split_\\d+/|metrics/|.csv', '', metric.files) %>% 
gsub('_metrics', '', .) %>%
gsub('/', '_', .)


metrics_df <- bind_rows(metrics, .id = 'celltype') %>%
mutate(
      gene.set = str_extract(celltype, "HVG\\.autosome|chrX|autosome|HVG|SLE"),
      split = str_extract(celltype, "split_[0-9]*"),
      age_group = str_extract(celltype, "adult|child"),
      celltype = gsub("^nehar_belaid_update_metrics_", "", celltype),
      celltype = gsub("_child_split_[0-9]*", "", celltype),
      celltype = gsub("\\.HVG\\.autosome|\\.SLE|\\.chrX|\\.autosome|\\.HVG", "", celltype),
      celltype = gsub("adult|child|split_[0-9]*", "", celltype),
      celltype = gsub("_", " ", celltype),
      celltype = gsub("  ", " ", celltype),
      celltype = gsub(" $", "", celltype)
    ) 
 
 
average_model_metrics_df <- metrics_df %>%
  group_by(celltype, gene.set, age_group) %>%
  summarise(
    Accuracy = mean(Accuracy),
    Precision = mean(Precision),
    Recall= mean(Recall),
    F1=mean(F1), 
    AUC=mean(AUC), 
    AUPRC=mean(AUPRC), 
    Kappa=mean(Kappa),
    MCC=mean(MCC), 
    n_features=round(mean(n_features),0)) %>% 
  data.frame()

 
 
   
  
  pdf("validation_mcc_combined_X_SLE_fig4.pdf")
  
  ggplot( metrics_df, aes(y=MCC, x=celltype, fill=gene.set)) +
    geom_boxplot(position=position_dodge(width=0.8), alpha=0.7) +
    geom_point(position=position_jitterdodge(dodge.width=0.8, jitter.width=0.2), 
               size=1) +
    geom_hline(yintercept = 0.7, linetype = 'dotted', color='red') + 
    geom_hline(yintercept = 0.5, linetype = 'dotted', color='grey') +
    facet_wrap(~age_group, ncol=1, nrow=5) + 
    theme_minimal() +
    labs(y='MCC', x='') +
    theme(
      axis.text.y=element_text(size=10),
      axis.text.x=element_text(size=14, angle=90),
      legend.position="top") +
    scale_fill_manual(values=gene.set.colours, name='Gene Set')
  dev.off() 
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  # qsub -I -l select=1:ncpus=2:mem=100gb -l walltime=02:00:00
  # cd /srv/scratch/ballouz/sle
  # module load r/4.4.0
  # R
  
  obj = readRDS("pbmc_female.control_managed.RDS")
  
  obj$cell_type_detailed2 = gsub("_|\\+|-", " ", obj$cell_type_detailed)
  load("colours.rdata") 
  obj$cell_type_detailed2 = factor(obj$cell_type_detailed2, levels = celltypes)
  
  
  metadata = obj@meta.data 
  save(metadata, file="metadata_pbmc_female.control_managed.Rdata")
  phens = unique( cbind(metadata$donor_id, metadata$age, metadata$sex,metadata$ancestry, metadata$disease) )
  save(phens, file="phens.Rdata")
  
  
  
  png("testumap.png", height=2000, width=2000)
  DimPlot(obj, 
          reduction = "umap", 
          pt.size = 2, 
          raster=FALSE, 
          group.by = "cell_type_detailed2" , 
          cols=  celltypes.colours ) + NoLegend() 
  dev.off()
  
  
  
  
  
  
  obj = readRDS("pbmc.female.RDS")
  
  
  load("colours.rdata") 
  celltypes_2 = unique(obj@meta.data$cellTypist)
  
  
  
  library(Seurat)
  library(dplyr)
  library(edgeR)
  
  # qsub -I -l select=1:ncpus=2:mem=100gb -l walltime=02:00:00
  
  
  Perez <- obj$cellTypist
  
  Perez <- case_when(
    Perez == "NK cells" ~ "natural killer cell",
    Perez == "Naive B cells" ~ "B cell Naive",
    Perez == "B cells" ~ "B cell",
    Perez == "Classical monocytes" ~ "Classical monocyte",
    Perez == "Non-classical monocytes" ~ "Non-classical monocyte",
    Perez == "Tcm/Naive helper T cells" ~ "CD4 T cell Naive",
    Perez == "Tcm/Naive cytotoxic T cells" ~ "CD8 T cell Naive",
    Perez == "Tem/Trm cytotoxic T cells" ~ "CD8 positive, alpha beta T cell",
    Perez == "Tem/Effector helper T cells" ~ "CD4 positive, alpha beta T cell",
    Perez == "Regulatory T cells" ~ "CD4 T cell Treg",
    Perez == "MAIT cells" ~ "CD8 T cell MAIT",
    Perez == "DC1" ~ "Conventional DC",
    Perez == "DC2" ~ "Conventional DC",
    Perez == "CD16+ NK cells" ~ "NK cell Bright",
    Perez == "pDC" ~ "Plasmacytoid DC",
    Perez == "Age-associated B cells" ~ "B cell Atypical",
    Perez == "Plasmablasts" ~ "Plasmablast",
    Perez == "HSC/MPP" ~ "Progenitor cell",
    Perez == "Tem/Temra cytotoxic T cells" ~ "CD8 T cell Cytotoxic GZMH ",
    #Perez == "Tem/Trm cytotoxic T cells" ~ "CD8 T cell Cytotoxic GZMK ",  
    TRUE ~ Perez
  )
  
  obj@meta.data$cell_type_detailed  <- Perez
  obj$cell_type_detailed2 = gsub("_|\\+|-", " ", obj$cell_type_detailed)
  
  
  metadata = obj@meta.data 
  save(metadata, file="metadata_pbmc.female.Rdata")
  phens = unique( cbind(metadata$donor_id, metadata$age, metadata$sex,metadata$ancestry, metadata$disease) )
  save(phens, file="phens.pbmc.female.Rdata")
  
  
  
  png("neharumap.png", height=2000, width=2000)
  DimPlot(obj, 
          reduction = "umap", 
          pt.size = 2, 
          raster=FALSE, 
          group.by = "cell_type_detailed2" , 
          cols=  celltypes.colours ) + NoLegend() 
  dev.off()
  
  
  names(celltypes.colours_expanded) =  c(names(celltypes.colours), setdiff(unique(obj$cell_type_detailed2), names(celltypes.colours)) )
  
  png("neharumap2.png", height=2000, width=2000)
  DimPlot(obj, 
          reduction = "umap", 
          pt.size = 2, 
          raster=FALSE, 
          group.by = "cell_type_detailed2" , 
          cols= celltypes.colours_expanded   ) + NoLegend() 
  dev.off()
  
  f1 = grepl("^c", metadata$condition)
  f2 = !f1
  
  
  png("neharumap_adult.png", height=2000, width=2000)
  DimPlot(obj[,f2], 
          reduction = "umap", 
          pt.size = 2, 
          raster=FALSE, 
          group.by = "cell_type_detailed2" , 
          cols= celltypes.colours_expanded   ) + NoLegend() 
  dev.off()
  
  
  
  png("neharumap_ped.png", height=2000, width=2000)
  DimPlot(obj[,f1], 
          reduction = "umap", 
          pt.size = 2, 
          raster=FALSE, 
          group.by = "cell_type_detailed2" , 
          cols= celltypes.colours_expanded   ) + NoLegend() 
  dev.off()
  
  
  
  library(ggplot2)
  library(reshape2)
  library(ggsignif)
  library(pwr)
  library(lsr)
  library(ROCR)
  
  
  
  
  df_orig <- read.delim('Flow_data/SLE_Lachie_Results.txt')
  df_orig$Disease <- factor(df_orig$Disease)
  df_orig$Run <- factor(df_orig$Run)
  df_orig$disease_activity <- ''
  df_orig$disease_activity[11:20] <- ifelse(df_orig$SLEDAI.2k.score[11:20] < 4, 'inactive', 'active')
  df_orig$disease_activity[1:10] <- 'healthy'
  df_orig$disease_activity <- factor(df_orig$disease_activity)
  
  df <- read.table('Flow_data/IL2RG_CD4_CD8.txt')
  df$Disease <- factor(df$CONDITION)
  df$AGE = df_orig$Age
  df$disease_activity <-  df_orig$disease_activity
  
  
  ggplot(df, aes(x=log.CD8_IL2RG., y=log.CD4_IL2RG., color=AGE) ) + geom_point() 
 
  comparisons <- list(c('SLE', 'HC'))
  
  
 pdf("il2rg_exp.pdf")  
  ggplot(df, aes(x=Disease, y= log.CD8_IL2RG.  )) +
    geom_boxplot(outlier.shape = NA,  fill= c('lightgreen', "yellow")) +
    geom_jitter(width = 0.2, size = 1.5, alpha = 0.6 ) +
    geom_signif(comparisons=comparisons, map_signif_level=TRUE, test='wilcox.test', color = 'black') +
    ylab('log(MFI)') + xlab('') +
    theme_minimal() + 
    theme(legend.position = "none")
  
  ggplot(df, aes(x=Disease, y= log.CD4_IL2RG.  )) +
    geom_boxplot(outlier.shape = NA,  fill= c('lightgreen', "yellow") ) +
    geom_jitter(width = 0.2, size = 1.5, alpha = 0.6) +
    geom_signif(comparisons=comparisons, map_signif_level=TRUE, test='wilcox.test', color = 'black') +
    ylab('log(MFI)') + xlab('') +
    theme_minimal() +
       theme(legend.position = "none")
  dev.off() 
  
  
  
  CD4_IL2RG_cor <- cor.test( df$log.CD4_IL2RG.[11:20], df$SLEDAI[11:20], method="s")
  CD8_IL2RG_cor <- cor.test( df$log.CD8_IL2RG.[11:20], df$SLEDAI[11:20], method="s")
  
  ggplot(df[11:20,], aes(y =log.CD8_IL2RG., x = SLEDAI)) +
    geom_point() +
    geom_smooth(method = "lm", se = FALSE) +
    labs(title = "IL2RG Expression in CD8+ vs SLEDAI-2K Score",
         subtitle = paste("Correlation:", round(CD8_IL2RG_cor$estimate, 2),
                          ';', "p-value:", round(CD8_IL2RG_cor$p.value, 3)),
         y = "log(MFI)", x = "SLEDAI-2K") +
    theme_minimal()
  
  
  ggplot(df[11:20,], aes(y =log.CD4_IL2RG., x = SLEDAI)) +
    geom_point() +
    geom_smooth(method = "lm", se = FALSE) +
    labs(title = "IL2RG Expression in CD4+ vs SLEDAI-2K Score",
         subtitle = paste("Correlation:", round(CD4_IL2RG_cor$estimate, 2),
                          ';', "p-value:", round(CD4_IL2RG_cor$p.value, 3)),
         y = "log(MFI)", x = "SLEDAI-2K") +
    theme_minimal()
  
  
  
  
  ggplot(df, aes(x=AGE, y= log.CD4_IL2RG. , colour=Disease )) +
     geom_point() + 
     ylab('log(MFI)') + xlab('') +
    theme_minimal() + 
    theme(legend.position = "none")
  
  
  
  ggplot(df, aes(y=log.CD8_IL2RG., x= AGE , colour=disease_activity )) +
    geom_point() + 
   
    theme_minimal() + 
    theme(legend.position = "top")
  
   
  