setwd('')

library(reshape2)
library(ggplot2)
library(ggVennDiagram)
library(ggrepel)
library(dplyr)
library(WGCNA)
library(colormap)
library(limma)
library(pheatmap)
library(tidyr)
library(enrichR)
library(Rsamtools)
library(GenomicFeatures)
library(GenomicAlignments)
library(BiocParallel)
#start from raw counts
counts_matrix = read.csv('full counts matrix CTRP13.csv')
row.names(counts_matrix) = counts_matrix$gene_symbol
counts_matrix$gene_symbol = NULL
counts_matrix$X=NULL
ENS_ID = read.delim('MGI_all_gene_annot_with ENSMB.txt')
counts_matrix$gene_symbol = ENS_ID$X3..marker.symbol[(match(row.names(counts_matrix), ENS_ID$X11..Ensembl.gene.id))]
counts_matrix = counts_matrix[!duplicated(counts_matrix$gene_symbol),]
counts_matrix = counts_matrix[!is.na(counts_matrix$gene_symbol),]

row.names(counts_matrix)  = counts_matrix$gene_symbol
counts_matrix$gene_symbol=NULL
#52,183 genes
cc1 = counts_matrix[rowSums(counts_matrix) > 50,]
#23,395 genes


melted_cnts = melt(as.matrix(counts_matrix))
head(melted_cnts)
colnames(melted_cnts) = c('gene_symbol', 'SampleID', 'counts')
melted_cnts$genotype = ifelse(grepl('WT', melted_cnts$SampleID), 'WT', 'KO')
melted_cnts$tissue = gsub('KO_', '', melted_cnts$SampleID)
melted_cnts$tissue = gsub('WT_', '', melted_cnts$tissue)
melted_cnts$tissue = gsub('.bam', '', melted_cnts$tissue, fixed = T)
melted_cnts$tissue = gsub('[0-9]+', '', melted_cnts$tissue)
melted_cnts$tissue = gsub('G', 'GWAT', melted_cnts$tissue)
melted_cnts$tissue = gsub('i', 'iWAT', melted_cnts$tissue)
melted_cnts$tissue = gsub('L', 'Liver', melted_cnts$tissue)
melted_cnts$tissue = gsub('M', 'Muscle', melted_cnts$tissue)
melted_cnts$animal_ID =  extract_numeric(melted_cnts$SampleID)
melted_cnts$animal_ID = paste0(melted_cnts$genotype, '_', melted_cnts$animal_ID)
melted_cnts$logcnts = log2(melted_cnts$counts + 1)
melted_cnts$gene_tissue = paste0(melted_cnts$gene_symbol, '_', melted_cnts$tissue)

new_cnts = dcast(melted_cnts, gene_tissue ~ animal_ID, value.var = 'logcnts', fun.aggregate = mean, na.rm=T)
row.names(new_cnts) = new_cnts$gene_tissue
new_cnts$gene_tissue=NULL
new_cnts[1:10,1:10]
cnts_mat = as.data.frame(t(new_cnts))

cnts_mat$dm = ifelse(grepl('WT', row.names(cnts_mat)), 'WT', 'KO')
cnts_mat$dm = factor(cnts_mat$dm, levels=c('WT', 'KO'))
table(cnts_mat$dm)
design = model.matrix(~dm, data=cnts_mat)
head(design)
table(cnts_mat$dm)
dim(design)
new_cnts1 = as.data.frame(t(cnts_mat[, !colnames(cnts_mat)=='dm']))
fit = lmFit(new_cnts1, design)
fit = eBayes(fit)
row.names(fit)[1:10]

res_table = topTable(fit, coef=NULL,number=Inf, genelist=row.names(fit), adjust.method="BH", sort.by="B", resort.by=NULL, p.value=1, lfc=0, confint=FALSE)


write.csv(res_table, file = 'results from limma on KO over WT.csv', row.names = F)

res1 = res_table[!res_table$AveExpr<0.6,]
head(res1)
#need to play around to assessing proper thresholds.  
test_genes = as.vector(res1$ID[res1$P.Value<0.00005])
label_key = res1$ID[res1$P.Value<0.00005]
res1$label2 = ifelse(res1$ID %in% label_key, paste0(res1$ID), '')
table(res1$label2)
color_key_table = as.data.frame(res1$ID)
colnames(color_key_table) = 'ID'
color_key_table$tissue =gsub(".*_", "", color_key_table$ID)
color_sets = as.data.frame(table(color_key_table$tissue))

colnames(color_sets) = c('tissue_ID', 'Freq')
color_sets$colors = c('deepskyblue2', 'darkorange2', 'darkorchid2', 'firebrick')
color_key_table$color = color_sets$colors[match(color_key_table$tissue, color_sets$tissue_ID)]

res1$label_col1 = color_key_table$color[match(res1$ID, color_key_table$ID)]
res1$label_col2 = ifelse(res1$P.Value<0.001, paste0(res1$label_col1), 'gray74')

res1$tissue = gsub(".*_", "", res1$ID)
res1$gene_symbol = gsub("\\_.*","",res1$ID)
#Number of genes which will be labelled
#Volcano plot
pdf(file = 'Volcano Plot of KO over WT.pdf')
ggplot(res1, aes(x=logFC, y=-log10(P.Value))) + theme_classic() +
  geom_point(aes(x=logFC, y=-log10(P.Value)), color=res1$label_col2) +
  geom_label_repel(aes(x=logFC, y=-log10(P.Value), label = res1$label2), color = res1$label_col2, size = 2, label.size=NA, box.padding = 0.8, point.padding = 0.5, max.overlaps = Inf, alpha = .6, segment.color = 'grey50')  +   ggtitle('Volcano plot KO over WT')
dev.off()

new_pal = color_sets$colors
names(new_pal) = color_sets$tissue_ID
pdf(file = paste0('TISSUE LEGEND for Volcano Plot of KO over WT.pdf'))
plot.new()
legend("center",
       legend = names(new_pal),
       fill = new_pal,       # Color of the squares
       border = "black")
dev.off()

proportions_plot = function(pval_threshold) {
  cof1_table = res1[res1$P.Value<pval_threshold,]
  cof1_table$FC_cat = ifelse(cof1_table$logFC>1.2, 'up-regulated', 'down-regulated')
  cof1_table$tissue = color_key_table$tissue[match(cof1_table$ID, color_key_table$ID)]
  
  pp1 = cof1_table %>%
    group_by(tissue) %>%
    summarise(cnt = n()) %>%
    mutate(freq = round(cnt / sum(cnt), 3)) %>% 
    arrange(desc(freq))
  gene_count = length(row.names(cof1_table))
  
  pdf(file = paste0('Proportions of KO gene changes P_less ', pval_threshold, '.pdf'))
  g1 = ggplot(pp1, aes(x=tissue, y=freq, fill=tissue)) + geom_bar(stat="identity", width=.5, position = "dodge") + theme_classic() + scale_fill_hue(l=20, c=60) + ylab('Frequency by datatype') + xlab('')+ theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1)) + ggtitle(paste0('KO changes P < ', pval_threshold,  ', ', gene_count, ' genes'))
  print(g1)
  dev.off()
}
proportions_plot(0.001)
  
  res1$tissue = gsub(".*_", "", res1$ID)
  res1$gene_symbol = gsub("\\_.*","",res1$ID)
  table(res1$tissue)
  #intersect genes

  intersect_degs = function(pval_threshold){
  cof1_table = res1[res1$P.Value<pval_threshold,]
  gene_list = list(GWAT = cof1_table$gene_symbol[cof1_table$tissue=='GWAT'], iWAT = cof1_table$gene_symbol[cof1_table$tissue=='iWAT'], Muscle = cof1_table$gene_symbol[cof1_table$tissue=='Muscle'], Liver = cof1_table$gene_symbol[cof1_table$tissue=='Liver'])
  pdf(file = paste0('Intersection of DEGs pvalue less ', pval_threshold, '.pdf'))
  g1 = ggVennDiagram(gene_list, label_alpha = 0) + scale_fill_distiller(palette = "PiYG") + ggtitle(paste0('Intersection of DEGs, pvalue < ', pval_threshold))
  print(g1)
  dev.off()
  }
  
intersect_degs(0.001)




#####################
#tissue_enrichments
setEnrichrSite("Enrichr")
dbs1 <- c("GO_Biological_Process_2021", "KEGG_2021_Human", "Reactome_2022")

tissue_run = 'iWAT'

tissie_enrichments = function(tissue_run){
  pp1 = res1[res1$tissue %in% tissue_run,]
pp1 = pp1[pp1$P.Value < 0.01,]
pp1 = pp1[pp1$logFC>1.2,]
pp1_length = ifelse(length(row.names(pp1)) > 500, as.numeric(500), as.numeric(length(row.names(pp1))))
pp2 = pp1[1:pp1_length,]
gg1 = pp2$gene_symbol

enriched <- enrichr(gg1, dbs1)
pdf(file = paste0('CTRP13 KO upregulated genes', tissue_run, ' ', names(enriched[1]), '.pdf'))
g1 = plotEnrich(enriched[[1]], showTerms = 15, numChar = 40, y = "Count", orderBy = "P.value") + ggtitle(paste0('CTRP13 KO upregulated genes', tissue_run, ' ', names(enriched[1])))
print(g1)
dev.off()

pdf(file = paste0('CTRP13 KO upregulated genes', tissue_run, ' ', names(enriched[2]), '.pdf'))
g2 = plotEnrich(enriched[[2]], showTerms = 15, numChar = 40, y = "Count", orderBy = "P.value") + ggtitle(paste0('CTRP13 KO upregulated genes', tissue_run, ' ', names(enriched[2])))
print(g2)
dev.off()

pdf(file = paste0('CTRP13 KO upregulated genes', tissue_run, ' ', names(enriched[3]), '.pdf'))
g3 = plotEnrich(enriched[[3]], showTerms = 15, numChar = 40, y = "Count", orderBy = "P.value") + ggtitle(paste0('CTRP13 KO upregulated genes', tissue_run, ' ', names(enriched[3])))
print(g3)
dev.off()

pp1 = res1[res1$tissue %in% tissue_run,]
pp1 = pp1[pp1$P.Value < 0.01,]
pp1 = pp1[pp1$logFC<1.2,]
pp1_length = ifelse(length(row.names(pp1)) > 500, as.numeric(500), as.numeric(length(row.names(pp1))))
pp2 = pp1[1:pp1_length,]
gg1 = pp2$gene_symbol

enriched <- enrichr(gg1, dbs1)
pdf(file = paste0('CTRP13 KO downregulated genes', tissue_run, ' ', names(enriched[1]), '.pdf'))
g1 = plotEnrich(enriched[[1]], showTerms = 15, numChar = 40, y = "Count", orderBy = "P.value") + ggtitle(paste0('CTRP13 KO downregulated genes', tissue_run, ' ', names(enriched[1])))
print(g1)
dev.off()

pdf(file = paste0('CTRP13 KO downregulated genes', tissue_run, ' ', names(enriched[2]), '.pdf'))
g2 = plotEnrich(enriched[[2]], showTerms = 15, numChar = 40, y = "Count", orderBy = "P.value") + ggtitle(paste0('CTRP13 KO downregulated genes', tissue_run, ' ', names(enriched[2])))
print(g2)
dev.off()

pdf(file = paste0('CTRP13 KO downregulated genes', tissue_run, ' ', names(enriched[3]), '.pdf'))
g3 = plotEnrich(enriched[[3]], showTerms = 15, numChar = 40, y = "Count", orderBy = "P.value") + ggtitle(paste0('CTRP13 KO downregulated genes', tissue_run, ' ', names(enriched[3])))
print(g3)
dev.off()
}
tissie_enrichments('GWAT')
tissie_enrichments('iWAT')
tissie_enrichments('Muscle')
tissie_enrichments('Liver')



###############################################################
#METSIM Analysis
#read in metsim data
metsim_genes = read.delim('METSIM_trx_adipose.txt')
metsim_clinical = read.delim('METSIM_clinical.txt')
sig_degs = res1[res1$P.Value<0.0005 & res1$tissue=='GWAT',]
orths = read.delim('raw_data/Mouse Gene info with Human Orthologues.txt')
sig_degs$human_orth = orths$human_orth[match(sig_degs$gene_symbol, orths$Symbol)]
ctrp13 = metsim_genes[metsim_genes$gene_symbol=='C1QL3',]
ctrp13 = na.omit(ctrp13)

desc_set = c('B: ApoB', 'B: Bioimpedance: Basal metabolic rate (kcal)', 'B: Bioimpedance: Fat mass (%)', ' B: IL1 beta (pg/ml)', 'B: LDL cholesterol (mmol/l)', 'B: OGTT fasting plasma FFA (mmol/l)', 'B: OGTT fasting plasma insulin (mU/l)', 'B: OGTT fasting plasma glucose (mmol/l)', 'B: Plasma glucose area under the curve (OGTT) (mmol/l * min)', 'B: Total triglycerides (mmol/l)', 'B: Waist/hip ratio')
gset = c('C1QL3', sig_degs$human_orth[sig_degs$logFC>1.2])
mm1 = metsim_genes[metsim_genes$gene_symbol %in% gset,]
mm2 = metsim_clinical[metsim_clinical$trait_description %in% desc_set,]
mm3 = dcast(mm1, metsim_id ~gene_symbol, value.var='expression_value', fun.aggregate=mean, na.rm=T)
mm4 = dcast(mm2, metsim_id ~trait_description, value.var='trait_value', fun.aggregate=mean, na.rm=T)
row.names(mm3) = mm3$metsim_id
mm3$metsim_id=NULL
row.names(mm4) = mm4$metsim_id
mm4$metsim_id=NULL
mm4 = mm4[row.names(mm4) %in% row.names(mm3),]
mm3 = mm3[row.names(mm3) %in% row.names(mm4),]

cc1 = bicorAndPvalue(mm3, mm4, use = 'p')
gset1 = melt(cc1$bicor)
gset1$pvalue = melt(cc1$p)$value
gset1 = gset1[gset1$pvalue<0.01 & gset1$value<0,]
gset1 = na.omit(gset1)
upreg_neg = as.vector(unique(gset1$Var1))

gset = c('C1QL3', sig_degs$human_orth[sig_degs$logFC<1.2])
mm1 = metsim_genes[metsim_genes$gene_symbol %in% gset,]
mm2 = metsim_clinical[metsim_clinical$trait_description %in% desc_set,]
mm3 = dcast(mm1, metsim_id ~gene_symbol, value.var='expression_value', fun.aggregate=mean, na.rm=T)
mm4 = dcast(mm2, metsim_id ~trait_description, value.var='trait_value', fun.aggregate=mean, na.rm=T)
row.names(mm3) = mm3$metsim_id
mm3$metsim_id=NULL
row.names(mm4) = mm4$metsim_id
mm4$metsim_id=NULL
mm4 = mm4[row.names(mm4) %in% row.names(mm3),]
mm3 = mm3[row.names(mm3) %in% row.names(mm4),]

cc1 = bicorAndPvalue(mm3, mm4, use = 'p')
gset1 = melt(cc1$bicor)
gset1$pvalue = melt(cc1$p)$value
gset1 = gset1[gset1$pvalue<0.01 & gset1$value<0,]
gset1 = na.omit(gset1)
downreg_pos = as.vector(unique(gset1$Var1))

full_cor_set = c('C1QL3', downreg_pos, upreg_neg)

mm1 = metsim_genes[metsim_genes$gene_symbol %in% full_cor_set,]
mm2 = metsim_clinical[metsim_clinical$trait_description %in% desc_set,]

mm3 = dcast(mm1, metsim_id ~gene_symbol, value.var='expression_value', fun.aggregate=mean, na.rm=T)
mm4 = dcast(mm2, metsim_id ~trait_description, value.var='trait_value', fun.aggregate=mean, na.rm=T)
row.names(mm3) = mm3$metsim_id
mm3$metsim_id=NULL
row.names(mm4) = mm4$metsim_id
mm4$metsim_id=NULL
mm3 = mm3[row.names(mm3) %in% row.names(mm4),]
mm4 = mm4[row.names(mm4) %in% row.names(mm3),]

cc1 = bicorAndPvalue(mm3, mm4, use = 'p')
cc2 = cc1$bicor
cc2[is.na(cc2)] = 0
tt4 = cc1$p
tt4[is.na(tt4)] = 1

#heatmap plot insert



####################
#Adjustment cors

sig_degs = res1[res1$P.Value<0.0001 & res1$tissue=='GWAT',]
orths = read.delim('G:/My Drive/Datasets/Mouse/genome files/Mouse Gene info with Human Orthologues.txt')
sig_degs$human_orth = orths$human_orth[match(sig_degs$gene_symbol, orths$Symbol)]
ctrp13 = metsim_genes[metsim_genes$gene_symbol=='C1QL3',]
ctrp13 = na.omit(ctrp13)


desc_set = c('B: ApoB', 'B: Bioimpedance: Basal metabolic rate (kcal)', 'B: Bioimpedance: Fat mass (%)', ' B: IL1 beta (pg/ml)', 'B: LDL cholesterol (mmol/l)', 'B: OGTT fasting plasma FFA (mmol/l)', 'B: OGTT fasting plasma insulin (mU/l)', 'B: OGTT fasting plasma glucose (mmol/l)', 'B: Plasma glucose area under the curve (OGTT) (mmol/l * min)', 'B: Total triglycerides (mmol/l)', 'B: Waist/hip ratio')
gset = c('C1QL3', sig_degs$human_orth)
mm1 = metsim_genes[metsim_genes$gene_symbol %in% gset,]
mm2 = metsim_clinical[metsim_clinical$trait_description %in% desc_set,]

mm3 = dcast(mm1, metsim_id ~gene_symbol, value.var='expression_value', fun.aggregate=mean, na.rm=T)
mm4 = dcast(mm2, metsim_id ~trait_description, value.var='trait_value', fun.aggregate=mean, na.rm=T)
row.names(mm3) = mm3$metsim_id
mm3$metsim_id=NULL
row.names(mm4) = mm4$metsim_id
mm4$metsim_id=NULL
mm3 = mm3[row.names(mm3) %in% row.names(mm4),]
mm4 = mm4[row.names(mm4) %in% row.names(mm3),]
cc1 = bicorAndPvalue(mm3$TNFRSF1B, mm4, use = 'p')

pstat_table = melt(cc1$p)
colnames(pstat_table) = c('gene_symbol', 'trait_name', 'pvalue')
pstat_table$category = paste0('un-adjusted')
resids = summary(lm(mm3$C1QL3 ~ mm3$TNFRSF1B))$residuals
resid_cors = bicorAndPvalue(resids, mm4, use = 'p')
adj_table = melt(resid_cors$p)
colnames(adj_table) = c('gene_symbol', 'trait_name', 'pvalue')
adj_table$category = paste0('C1QL3-adjusted')

full_table = as.data.frame(rbind(pstat_table, adj_table))
full_table$logp = -log10(full_table$pvalue)
full_table$category = factor(full_table$category, levels = c('un-adjusted', 'C1QL3-adjusted'))
pdf(file = 'TNFRSF1B correlations with traits with C1QL3 mediation.pdf')
ggplot(full_table, aes(x=trait_name, y=logp, fill = category)) + geom_col(position = "dodge") + theme_minimal() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + ylab('-log10(pvalue_ of regression') + scale_fill_manual(values = c('dodgerblue2', 'darkorange3')) + xlab('') + ggtitle('TNFRSF1B correlations with traits with C1QL3 mediation')
dev.off()



cc1 = bicorAndPvalue(mm3$RELN, mm4, use = 'p')

pstat_table = melt(cc1$p)
colnames(pstat_table) = c('gene_symbol', 'trait_name', 'pvalue')
pstat_table$category = paste0('un-adjusted')
resids = summary(lm(mm3$C1QL3 ~ mm3$RELN))$residuals
resid_cors = bicorAndPvalue(resids, mm4, use = 'p')
adj_table = melt(resid_cors$p)
colnames(adj_table) = c('gene_symbol', 'trait_name', 'pvalue')
adj_table$category = paste0('C1QL3-adjusted')

full_table = as.data.frame(rbind(pstat_table, adj_table))
full_table$logp = -log10(full_table$pvalue)
full_table$category = factor(full_table$category, levels = c('un-adjusted', 'C1QL3-adjusted'))
pdf(file = 'RELN correlations with traits with C1QL3 mediation.pdf')
ggplot(full_table, aes(x=trait_name, y=logp, fill = category)) + geom_col(position = "dodge") + theme_minimal() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + ylab('-log10(pvalue_ of regression') + scale_fill_manual(values = c('dodgerblue2', 'darkorange3')) + xlab('') + ggtitle('RELN correlations with traits with C1QL3 mediation')
dev.off()
