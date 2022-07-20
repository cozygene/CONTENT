############################### Hierarchical multiple testing for CONTENT
############################### Brunilda Balliu and Mike Thompson
############################### Los Angeles, November 17th, 2020

# The get_eGenes_multi_tissue command below is going to save a bunch of eAssociation files to the working directory (!!!). 
# Set the place you want them to be saved by changing the working directory 
# setwd(dir = 'GTEx_v7_rsq_pval_20201117/') 
# Then, make a folder called "m_eqtl_out_tmp" inside that folder to save some temporary files (see below)

############# Libraries
library(data.table) # make sure to have the newest version of this package
library(TreeQTL) #make sure to download the newest version from http://www.bioinformatics.org/treeqtl/
library(dplyr)
library(ggplot2)
setwd("~/content_paper/finalfdr/content_ldref/") # change this to where your summary statistics are
############# R-square p-values from expression imputation
# each of these must follow the format specified in the zenodo link
# gene tissue pvalue r2 tab-delimited columns
full=fread(file = "content_full_ldref.txt",header = T)
specific=fread(file = "content_specific_ldref.txt",header = T)
shared=fread(file = "content_shared_ldref.txt",header = T) %>% mutate(tissue="Average")

############# Arguments
# tissues =  Average, Adipose_Subcutaneous, Adipose_Visceral_Omentum, etc
# snps = Model, i.e. Full, shared, specific
genes_names=union(union(full$gene,specific$gene),shared$gene) 
tissue_names = sort(union(union(full$tissue,specific$tissue),shared$tissue)) #vector of names for each tissue in ****alphanumeric order***, "Average" = shared component
snp_names=c("Full","shared","specific")

# process files to create genes_by_tissue and snps_by_tissue data frames
full_wide=data.frame(reshape(data = full %>% select(-rsq), idvar = "gene", timevar = "tissue", direction = "wide"),row.names = 1,check.names = F)
colnames(full_wide) = gsub(pattern = "pvalue.",replacement = "",x = colnames(full_wide))
full_wide[!is.na(full_wide)]=1
full_wide[is.na(full_wide)]=0

shared_wide=data.frame(reshape(data = shared %>% select(-rsq), idvar = "gene", timevar = "tissue", direction = "wide"),row.names = 1,check.names = F) 
colnames(shared_wide) = gsub(pattern = "pvalue.",replacement = "",x = colnames(shared_wide))
shared_wide[!is.na(shared_wide)]=1
shared_wide[is.na(shared_wide)]=0

genes_by_tissue=merge(x = data.table(full_wide,keep.rownames = T),y = data.table(shared_wide,keep.rownames = T)) %>% setnames(old = "rn", new = "genes") %>% select(c("genes",all_of(tissue_names))) %>% data.frame(check.names = F)

snps_by_tissue=data.frame(snp=snp_names, matrix(data = 1,nrow = 3, ncol = length(tissue_names),dimnames = list(NULL,tissue_names)),check.names = F)

gene_map = data.frame(gene = genes_names, chromosome=1, start=1, end=2)
snp_map = data.frame(snp = snp_names, chromosome=1, bp=1)

m_eqtl_out_dir = "m_eqtl_out_tmp/" # temporary directory where MatrixEQTL outputs will be saved

alpha_level1 = .05
alpha_level2 = .05
alpha_level3 = .05


############# Create temporary MatrixEQTL outputs
for(tissue_i in tissue_names){
  multi_tissue_pvals = rbind(
    full %>% filter(tissue==tissue_i) %>%  mutate(SNP="Full", beta=NA, `t-stat`=NA,  FDR = NA)  %>% setnames(old = "pvalue", new = "p-value") %>% select(SNP, gene, beta, `t-stat`, `p-value`,  FDR) %>% filter(!is.na(`p-value`)), 
    
    shared %>% filter(tissue==tissue_i) %>%  mutate(SNP="shared", beta=NA, `t-stat`=NA,  FDR = NA)  %>% setnames(old = "pvalue", new = "p-value") %>% select(SNP, gene, beta, `t-stat`, `p-value`,  FDR) %>% filter(!is.na(`p-value`)),
    
    specific %>% filter(tissue==tissue_i) %>%  mutate(SNP="specific", beta=NA, `t-stat`=NA,  FDR = NA)  %>% setnames(old = "pvalue", new = "p-value") %>% select(SNP, gene, beta, `t-stat`, `p-value`,  FDR) %>% filter(!is.na(`p-value`)))
  
  
  # You can threshold the p-values that are very large if it takes too long to write and read these files. You do not need them to perform the correction.
  # multi_tissue_pvals <- multi_tissue_pvals[multi_tissue_pvals$`p-value` <= 0.1, ] 
  
  # Sort appropriately, not needed but just to match MatrixEQTL output format
  multi_tissue_pvals <- multi_tissue_pvals[order(multi_tissue_pvals$`p-value`), ]
  
  # Write to disk
  write.table(x = multi_tissue_pvals, file = paste0(m_eqtl_out_dir,tissue_i,"_rsq_pvals.txt"), append = F, sep = "\t", quote = FALSE, row.names = FALSE)
}

#############  Identify eGenes, i.e. heritable genes in at least one tissue according to at least one model
# Takes a couple of minutes to run, saves files called eAssoc_3level_by_gene_tissue_TISSUE_NAME.txt inside the working directory
eGenes=get_eGenes_multi_tissue(genes_by_tissue = genes_by_tissue,snps_by_tissue = snps_by_tissue, 
                               gene_map = gene_map, snp_map = snp_map,
                               nearby = T,dist = 1e06,
                               m_eqtl_out_dir = m_eqtl_out_dir, tissue_names = tissue_names,
                               level1 = alpha_level1, level2 = alpha_level2 ,level3 = alpha_level3)
write.table(x = eGenes,file = "eGenes_3level_by_gene_tissue.txt", append = F, sep = "\t", quote = FALSE, row.names = FALSE)

# Clean up
unlink(paste0(m_eqtl_out_dir,"*","_rsq_pvals.txt"))

############# Merge files together 
all_eAssoc=NULL

for(tissue_i in tissue_names){
  eAssoc_tissue_i = fread(input = paste0("eAssoc_3level_by_gene_tissue_",tissue_i,".txt")) %>% mutate(tissue=tissue_i)
  all_eAssoc=rbind(all_eAssoc,eAssoc_tissue_i)
}
saveRDS(all_eAssoc,"contentassoc_ldref.rds")
head(all_eAssoc)

