#!/usr/bin/env Rscript
suppressPackageStartupMessages(require(optparse))

option_list = list(
  make_option(c("-X", "--X_file"), action="store", default=NULL, type='character',
              help="A genotype matrix with the set of individuals in all Y_files and their corresponding
              cis-SNPs. Individuals are row names, no column names."),
  make_option(c("-Y", "--Y_file_dir"), action="store", default=NULL, type='character',
              help="Directory containing only the expression files for each context. The individuals
              (row names) must be a subset of those in the genotype (X_file). Contains an unnamed
              column with the expression (or if cov_file_dir is NULL, the residual expression
              for a given context."),
  make_option(c("-c", "--cov_file_dir"), action="store", default=NULL, type='character',
              help="Specify this directory if you want to residualize the expression. The
              covariates are context-specific and should share the same prefix as the expression 
              files such that file i of the expression corresponds to file i of the covariates,
              e.g. Liver_expression and Liver_covariates. If NULL, CONTENT will run assuming
              the expression has already been residualized. Needs to be the same row-length as the
              corresponding expression file, contain the same rownames, and contain no column names."),
  make_option(c("-v", "--verbose"), action="store_true", dest="verbose", default=T,
              help="Make the program be verbose (show some progress/output.)"),
  make_option(c("-o", "--out_dir"), action="store", default=NULL, type='character',
              help="Where to save the cross-validated predictors and their corresponding performance
              statistics. Use the same directory for all genes, and it will contain the predictors as
              well as pvalues and rsquares for each gene for convenient reporting. E.g. gene1_pvals,
              gene1rsqs, gene2pvals, gene2rsqs"),
  make_option(c("-s", "--snps"), action="store", default=NULL, type='character',
               help="A tab-delimited txt file containing information from your .bed or other genotype
               file. Contains 6 columns and the number of rows corresponds to the number of snps. 
               No column names or row names. Col2 must be rsIDs for TWAS to work. 
               Col1 chromosome
               Col2 rsID
               Col3 location CM (this doesn't really matter for TWAS)
               Col4 location/locus on chromosome
               Col5 allele1
               Col6 allele2"),
  make_option(c("-t", "--twas_dir"), action="store", default=NULL, type='character',
               help="Where to save the TWAS weights. Two directories will spawn inside this directory,
               significant and nonsignificant corresponding to genes that were significantly predicted
               or not. Use the same dir for all genes and it will conveniently store the results (all
               significantly predicted genes will have their weights in the same directory so we can 
               just point TWAS to this single directory)."),
  make_option(c("-a", "--signif"), action="store", default=0.1, type='numeric',
               help="Nominally significant pvalue threshold (alpha). If a method/predictor significantly
               explains some variability of the gene expression at this threshold, save a TWAS weight 
               for this in the twas_dir/significant directory. Significant weights will be stored in
               twas_dir/significant and nonsignificant weights will be stored in twas_dir/not_significant
               in case you change your mind later on (check this using the saved performance statistics)."),
  make_option(c("-n", "--gene_name"), action="store", default=NULL, type='character',
               help="Add this prefix to all the saved results files... necessary to distinguish results
               where more than one gene-analysis is run."),
  make_option(c("-l", "--content_alpha"), action="store", default=0.5, type='numeric',
              help="The regularization constant for CONTENT. Default is 1 (LASSO). Minimum value
              of 1e-4."),
  make_option(c("-e", "--tissue_alpha"), action="store", default=0.5, type='numeric',
              help="The regularization constant for the tissue by tissue approach. Default is .5 (eNet).
              Minimum value of 1e-4."),
  make_option(c("-f", "--num_folds"), action="store", default=10, type='numeric',
              help="Number of folds for cross-validation."),
  make_option(c("-d", "--scale_exp"), action="store", dest="scale_exp", default=F,
              help="Center and scale the gene expression.")
)
opt = parse_args(OptionParser(option_list=option_list))

library(bigstatsr)
library(data.table)
library(caret)

X=NULL
Ys=NULL
covs=NULL
X_file=opt$X_file
Y_file_dir=opt$Y_file_dir
cov_file_dir=opt$cov_file_dir
verbose=opt$verbose
out_dir=opt$out_dir
snps=opt$snps
twas_dir=opt$twas_dir
gene_name=opt$gene_name
signif=opt$signif
num_folds=10
content_alpha=opt$content_alpha
tissue_alpha=opt$tissue_alpha
scale_exp=opt$scale_exp
seed=9000
## Open data and get directories set up
set.seed(seed)
if(verbose){
  message("Saving cross-validated predictors and performance metrics in ", out_dir)
}
if(is.null(twas_dir)){
  twas_dir=paste0(out_dir,"TWAS_weights")
  if(!dir.exists(twas_dir)){
    system(paste0("mkdir ", twas_dir)) 
  }
  if(verbose){
    message("Saving TWAS weights in ", twas_dir) 
  }
}else{
  if(verbose){
    message("Saving TWAS weights in ", twas_dir) 
  }
}
## By default, try to run content, but if there is only one context, do not
runcontent=T
# Check if X is a matrix or a path
## If matrix, just make sure it's good to go
if(class(X) %in% c("matrix", "data.frame")){
  if(verbose){
    message("Converting objects to matrix...")
  }
  X<-as.matrix(X)
  Ys<-lapply(Ys, as.matrix)
  if(length(Ys)<3){
    runcontent=F
  }
  covs<-lapply(covs, as.matrix)
  ## If X is a path, then everything is else treated as such
}else if(class(X_file) == "character"){
  if(verbose){
    message("Reading in files...") 
  }
  suppressWarnings(expr = {X<-fread(file = X_file, sep='\t', data.table=F)})
  X<-as.matrix(data.frame(X, row.names=1, check.names = F))
  if(!is.null(cov_file_dir)){
    covs<-vector("list", length = length(list.files(cov_file_dir)))
    names(covs)=list.files(cov_file_dir)
    for(i in 1:length(covs)){
      suppressWarnings(expr = {covs[[i]]<-fread(file = paste0(cov_file_dir,list.files(cov_file_dir)[i]), 
                                                sep='\t', data.table=F)})
      dts=sapply(covs[[i]], class)
      if (any(dts == "character") & min(which(dts == "character") != 1)) {
        stop("Covariate(s) ",
             paste0(which(dts == "character"), collapse=", "),
             " is(are) character(s). Please convert to numeric using one-hot
             encoding or as an integer.")
      }
      covs[[i]]<-as.matrix(data.frame(covs[[i]], row.names=1, check.names = F))
      }
  }
  Ys<-vector("list", length = length(list.files(Y_file_dir)))
  if(length(Ys)<3){
    runcontent=F
    message("Only 1 context found. CONTENT will not be run.")
  }
  names(Ys)<-list.files(Y_file_dir)
  for(i in 1:length(Ys)){
    suppressWarnings(expr={ Ys[[i]]<-fread(file = 
                                             paste0(Y_file_dir,list.files(Y_file_dir)[i]), sep='\t', data.table=F)})
    Ys[[i]]<-as.matrix(data.frame(Ys[[i]], row.names=1, check.names = F))
    if(any(is.na(Ys[[i]])) | any(is.nan(Ys[[i]]))){
      remove=unique(c( which(is.na(Ys[[i]])), which(is.nan(Ys[[i]])) ))
      Ys[[i]]=Ys[[i]][-remove,,drop=F]
    }
    if(scale_exp){
      Ys[[i]]=scale(Ys[[i]])
    }
  }
  }else{
    message("X is neither a matrix nor a path. Please supply X as a matrix, dataframe, 
            or a path.")
  }
if(is.null(snps)){
  message("snps file is null, TWAS weights will need to be updated after CONTENT is run...")
}else if(class(snps) == "character"){
  snps=read.table(snps, sep='\t', check.names = F, stringsAsFactors = F)
  storesnps=snps
}
if(is.null(gene_name)){
  message("gene_name not supplied, using 'gene' as gene name.
          results will not be easily distinguishable.")
}

# define some commonly-used variables
q<-length(Ys)
m<-ncol(X)
N<-nrow(X)
Yhats_tiss<-vector("list",q)
Yhats_hom<-vector("list",q)
Yhats_het<-vector("list", q)
for(i in 1:q){
  Yhats_het[[i]]<-matrix(NA, ncol=1, nrow=nrow(Ys[[i]]), 
                         dimnames = list(rownames(Ys[[i]]), "pred"))
  Yhats_hom[[i]]<-matrix(NA, ncol=1, nrow=nrow(Ys[[i]]), 
                         dimnames = list(rownames(Ys[[i]]), "pred"))
  Yhats_tiss[[i]]<-matrix(NA, ncol=1, nrow=nrow(Ys[[i]]), 
                          dimnames = list(rownames(Ys[[i]]), "pred"))
}

if(is.null(covs) & is.null(cov_file_dir)){
  if(verbose){
    message("Continuing without covariates")
  }
  ## The Ys are already residualized
  hom_expr_mat<-matrix(NA, nrow = nrow(X), ncol=q)
  rownames(hom_expr_mat)<-rownames(X)
  for(i in 1:q){
    hom_expr_mat[rownames(Ys[[i]]),i]<-Ys[[i]]
  }
}else{
  if(verbose){
    message("Residualizing over covariates.")
  }
  ## If there are covariates, generate residuals:
  hom_expr_mat<-matrix(NA, nrow = nrow(X), ncol=q)
  rownames(hom_expr_mat)<-rownames(X)
  for(i in 1:q){
    #### Regress out covariates for each tissue,
    ## If user supplies one directory for covariates, make sure to match
    covidx=match(names(Ys)[i], names(covs))
    tmp<-lm(formula = Ys[[i]] ~ covs[[covidx]][rownames(Ys[[i]]), ])
    ## check the names of residuals, if there's not a perfect overlap
    ### with rownames(Ys[[i]]) they had missing covs and we remove
    if(all(rownames(Ys[[i]]) %in% names(tmp$residuals))){
      hom_expr_mat[rownames(Ys[[i]]),i]<-tmp$residuals
    }else{
      use_inds=rownames(Ys[[i]])[which(rownames(Ys[[i]]) %in% names(tmp$residuals))]
      Ys[[i]]=Ys[[i]][use_inds,,drop=F]
      hom_expr_mat[rownames(Ys[[i]]),i]<-tmp$residuals[rownames(Ys[[i]])]
    }
    
  }
}

all_missing<-names(rowMeans(hom_expr_mat, na.rm = T)[which(is.nan(rowMeans(hom_expr_mat, na.rm = T)))])
remove_inds<-which(rownames(hom_expr_mat) %in% all_missing)
# These people are not in any tissues, we can just remove them from the data
if(length(remove_inds)>0){
  if(verbose){
    message("These individuals have NA in all Ys or covs: ")
    message(paste0(all_missing, collapse = ","))
  }
  hom_expr_mat<-hom_expr_mat[-remove_inds,,drop=F]
  X<-X[-remove_inds,]
}

# prepare for 10-fold cross-validation
test_inds_idx<-caret::createFolds(y=rowMeans(hom_expr_mat,na.rm=T), k = num_folds)
if(length(test_inds_idx)<num_folds){
  while(length(test_inds_idx)<num_folds){
    test_inds_idx<-caret::createFolds(y=rowMeans(hom_expr_mat,na.rm=T), k = num_folds)
  }
}
test_inds_ids<-lapply(test_inds_idx, function(x) {rownames(X)[x]})
test_inds<-list()
for(fold in 1:num_folds){
  for(i in 1:q){
    if(fold ==1){
      test_inds[[i]]<-list()
    }
    test_inds[[i]][[fold]]<-which(rownames(Ys[[i]]) %in% test_inds_ids[[fold]])
  }
}
if(verbose){
  message("Starting cross-validation")
  message("CONTENT temporary file is ", paste0(out_dir,gene_name, "_content_tmp.bk"))
}
if(file.exists(paste0(out_dir,gene_name, "_content_tmp.bk"))){
  system(paste0("rm ", paste0(out_dir,gene_name, "_content_tmp.bk")))
}
explanatory=as_FBM(X, backingfile=paste0(out_dir,gene_name, "_content_tmp"))
# start cross-validation
for(cur_fold in 1:num_folds){
  
  train_inds_id<-setdiff(rownames(X), test_inds_ids[[cur_fold]])
  train_inds_id<-intersect(train_inds_id, rownames(hom_expr_mat))
  fold_hom_expr_mat<-hom_expr_mat[train_inds_id,,drop=F]
  
  ## Fit the homogeneous component
  if(runcontent){
    hom_fit<-big_spLinReg(X = explanatory, match(rownames(fold_hom_expr_mat), rownames(X)),
                          y.train=rowMeans(x=fold_hom_expr_mat, na.rm = T),
                          alphas=c(content_alpha),K=10, warn=F)
    
    hom_beta_vals<-unlist(summary(hom_fit)$beta[
      which.min(summary(hom_fit)$validation_loss)])[1:ncol(explanatory)]
    hom_int<-unlist(summary(hom_fit)$intercept[
      which.min(summary(hom_fit)$validation_loss)])
    if(sum(is.na(hom_beta_vals))>0){
      hom_beta_vals[which(is.na(hom_beta_vals))]<-0
    }
    if(length(hom_beta_vals)<dim(explanatory)[2]){
      new_betas<-rep(0, dim(explanatory)[2])
      new_betas[attr(hom_fit, "ind.col")]<-hom_beta_vals
      hom_beta_vals<-new_betas
    }
    
    ## Fit the heterogeneous components
    het_tiss_betas<-vector("list", q)
    het_tiss_ints<-vector("list", q)
  }
  
  ## Also fit a tissue by tissue approach in case
  ## CONTENT cannot generate a model
  tiss_betas<-vector("list", q)
  tiss_ints<-vector("list", q)
  for(j in 1:q){
    cur_train_inds<-setdiff(rownames(Ys[[j]]), test_inds_ids[[cur_fold]])
    cur_train_inds<-intersect(cur_train_inds, rownames(fold_hom_expr_mat))
    if(runcontent){
      ### Find out which individuals were not present in other tissues
      nan_names<-names(rowMeans(fold_hom_expr_mat[cur_train_inds,-j],
                                na.rm = T)[which(is.nan(rowMeans(fold_hom_expr_mat[cur_train_inds,-j], na.rm = T)))])
      ### It may happen that there are individuals in this training set who
      ### do not appear in any other tissues, or they only appear in the test
      ### set of other tissues. If that's the case, remove them
      if(length(nan_names) == 0){
        cur_hom_expr_mat<-fold_hom_expr_mat[cur_train_inds,,drop=F]
      }else{
        subset_inds<-cur_train_inds[!(cur_train_inds %in% nan_names)]
        cur_hom_expr_mat<-fold_hom_expr_mat[subset_inds,,drop=F]
      }
      ### First, keep only the train individuals that are in this tissue+other tissues
      
      if(nrow(cur_hom_expr_mat) < 15){
        het_tiss_betas[[j]]<-rep(NA, m)
        het_tiss_ints[[j]]<-NA # 0
        tiss_inds=rownames(fold_hom_expr_mat)[which(!is.na(fold_hom_expr_mat[,j]))]
        tiss_response=fold_hom_expr_mat[tiss_inds,j]
        if(length(tiss_response) < 15){
          tiss_betas[[j]]<- rep(NA, m)
          tiss_ints[[j]]<-NA 
          next
        }
        tiss_fit<-big_spLinReg(X = explanatory, ind.train = match(tiss_inds, rownames(X)),
                               y.train = tiss_response,K=10, alphas = c(tissue_alpha),warn=F)
        
        tiss_beta_vals<-unlist(summary(tiss_fit)$beta[
          which.min(summary(tiss_fit)$validation_loss)])[1:m]
        tiss_tiss_int<-unlist(summary(tiss_fit)$intercept[
          which.min(summary(tiss_fit)$validation_loss)])
        idx_remove<-c()
        if(any(is.na(tiss_beta_vals))){
          idx_remove<-which(is.na(tiss_beta_vals))
          tiss_beta_vals[idx_remove]<-0
        }
        if(length(tiss_beta_vals)<dim(explanatory)[2]){
          new_betas<-rep(0, dim(explanatory)[2])
          new_betas[attr(tiss_fit, "ind.col")]<-tiss_beta_vals
          tiss_beta_vals<-new_betas
        }
        tiss_betas[[j]]<-tiss_beta_vals
        tiss_ints[[j]]<-tiss_tiss_int
        next
      }
      
      
      
      ### tissue_j expression - mean(all tissues except j expression)
      het_response<-cur_hom_expr_mat[, j]-rowMeans(x=cur_hom_expr_mat,  na.rm = T)
      
      if(length(het_response) < 15){
        het_tiss_betas[[j]]<- rep(NA, m)
        het_tiss_ints[[j]]<-NA #0
        tiss_inds=rownames(fold_hom_expr_mat)[which(!is.na(fold_hom_expr_mat[,j]))]
        tiss_response=fold_hom_expr_mat[tiss_inds,j]
        if(length(tiss_response) < 15){
          tiss_betas[[j]]<- rep(NA, m)
          tiss_ints[[j]]<-NA 
          next
        }
        tiss_fit<-big_spLinReg(X = explanatory, ind.train = match(tiss_inds, rownames(X)),
                               y.train = tiss_response,K=10, alphas = c(tissue_alpha),warn=F)
        
        tiss_beta_vals<-unlist(summary(tiss_fit)$beta[
          which.min(summary(tiss_fit)$validation_loss)])[1:m]
        tiss_tiss_int<-unlist(summary(tiss_fit)$intercept[
          which.min(summary(tiss_fit)$validation_loss)])
        idx_remove<-c()
        if(any(is.na(tiss_beta_vals))){
          idx_remove<-which(is.na(tiss_beta_vals))
          tiss_beta_vals[idx_remove]<-0
        }
        if(length(tiss_beta_vals)<dim(explanatory)[2]){
          new_betas<-rep(0, dim(explanatory)[2])
          new_betas[attr(tiss_fit, "ind.col")]<-tiss_beta_vals
          tiss_beta_vals<-new_betas
        }
        tiss_betas[[j]]<-tiss_beta_vals
        tiss_ints[[j]]<-tiss_tiss_int
        next
      }
      
      het_fit<-big_spLinReg(X = explanatory, ind.train = match(rownames(cur_hom_expr_mat), rownames(X)),
                            y.train = het_response,K=10, alphas = c(content_alpha),warn=F)
      
      het_beta_vals<-unlist(summary(het_fit)$beta[
        which.min(summary(het_fit)$validation_loss)])[1:m]
      het_tiss_int<-unlist(summary(het_fit)$intercept[
        which.min(summary(het_fit)$validation_loss)])
      idx_remove<-c()
      if(any(is.na(het_beta_vals))){
        idx_remove<-which(is.na(het_beta_vals))
        het_beta_vals[idx_remove]<-0
      }
      if(length(het_beta_vals)<dim(explanatory)[2]){
        new_betas<-rep(0, dim(explanatory)[2])
        new_betas[attr(het_fit, "ind.col")]<-het_beta_vals
        het_beta_vals<-new_betas
      }
      het_tiss_betas[[j]]<-het_beta_vals
      het_tiss_ints[[j]]<-het_tiss_int
      
    }
    
    
    tiss_inds=rownames(fold_hom_expr_mat)[which(!is.na(fold_hom_expr_mat[,j]))]
    tiss_response=fold_hom_expr_mat[tiss_inds,j]
    # tiss_explanatory=X[tiss_inds,]
    # 
    # tiss_explanatory<-as_FBM(tiss_explanatory, 
    #                          backingfile = paste0(out_dir,
    #                                               "content_tmp_",cur_fold))
    
    tiss_fit<-big_spLinReg(X = explanatory, ind.train = match(tiss_inds, rownames(X)),
                           y.train = tiss_response,K=10, alphas = c(tissue_alpha),warn=F)
    
    tiss_beta_vals<-unlist(summary(tiss_fit)$beta[
      which.min(summary(tiss_fit)$validation_loss)])[1:m]
    tiss_tiss_int<-unlist(summary(tiss_fit)$intercept[
      which.min(summary(tiss_fit)$validation_loss)])
    idx_remove<-c()
    if(any(is.na(tiss_beta_vals))){
      idx_remove<-which(is.na(tiss_beta_vals))
      tiss_beta_vals[idx_remove]<-0
    }
    if(length(tiss_beta_vals)<dim(explanatory)[2]){
      new_betas<-rep(0, dim(explanatory)[2])
      new_betas[attr(tiss_fit, "ind.col")]<-tiss_beta_vals
      tiss_beta_vals<-new_betas
    }
    tiss_betas[[j]]<-tiss_beta_vals
    tiss_ints[[j]]<-tiss_tiss_int
    
  }
  for(i in 1:q){
    safe_test_inds<-intersect(rownames(Ys[[i]])[test_inds[[i]][[cur_fold]]], rownames(Ys[[i]]))
    if(length(safe_test_inds) < 1){
      next
    }
    if(runcontent){
      Yhats_hom[[i]][rownames(Ys[[i]])[test_inds[[i]][[cur_fold]]],]<-X[rownames(Ys[[i]])[test_inds[[i]][[cur_fold]]], ] %*% hom_beta_vals + hom_int
      Yhats_het[[i]][rownames(Ys[[i]])[test_inds[[i]][[cur_fold]]],]<-X[rownames(Ys[[i]])[test_inds[[i]][[cur_fold]]], ] %*% het_tiss_betas[[i]] + het_tiss_ints[[i]]
    }
    Yhats_tiss[[i]][rownames(Ys[[i]])[test_inds[[i]][[cur_fold]]],]<-X[rownames(Ys[[i]])[test_inds[[i]][[cur_fold]]], ] %*% tiss_betas[[i]] + tiss_ints[[i]]
  }
  if(verbose & cur_fold < num_folds){
    message("Finished fold: ", cur_fold, " of ", num_folds)
  }
}
names(Yhats_het)<-names(Ys)
names(Yhats_hom)<-names(Ys)
names(Yhats_tiss)<-names(Ys)
save(Yhats_het, Yhats_hom, Yhats_tiss, file = paste0(out_dir,gene_name,"_crossval_predictors"))

## Score each method
het_cv_pvals<-vector("list", q); het_cv_r2s<-vector("list", q)
hom_cv_pvals<-vector("list", q); hom_cv_r2s<-vector("list", q)
het_cv_pvals.herit<-vector("list", q); het_cv_r2s.herit<-vector("list", q)
hom_cv_pvals.herit<-vector("list", q); hom_cv_r2s.herit<-vector("list", q)
full_cv_pvals<-vector("list", q); full_cv_r2s<-vector("list", q)
tiss_cv_pvals<-vector("list", q); tiss_cv_r2s<-vector("list", q)
full_weights<-vector("list", q); het_scales=vector("list", q)

## first is the hom term heritable:
if(runcontent){
  hom.tmp=do.call(rbind, Yhats_hom)
  hom.tmp=hom.tmp[rownames(hom_expr_mat),,drop=F]
  baseline=lm(rowMeans(hom_expr_mat[rownames(hom.tmp),], na.rm=T) ~ 1)
  homfit=lm(rowMeans(hom_expr_mat[rownames(hom.tmp),], na.rm=T) ~ hom.tmp)
  hom_test_stat.herit<--2*(logLik(baseline))+2*logLik(homfit)
  hom_cv_pvals.herit=pchisq(hom_test_stat.herit,1,lower.tail=F)
  hom_cv_r2s.herit=summary(homfit)$adj.r.squared
}else{
  hom_cv_pvals.herit=NA
  hom_cv_r2s.herit=NA
}

for(i in 1:q){
  
  # baseline model
  m1<-lm(hom_expr_mat[rownames(Ys[[i]]),i]~1)
  # homogeneous model
  if(runcontent){
    m2<-lm(hom_expr_mat[rownames(Ys[[i]]),i] ~ Yhats_hom[[i]][rownames(Ys[[i]]),])
  }
  # heterogeneous model
  ## if we learned one for this tissue/context
  if(!all(is.na(Yhats_het[[i]]))){
    ## is the het term heritable:
    hetresponse=hom_expr_mat[rownames(Ys[[i]]),i] - rowMeans(hom_expr_mat[rownames(Ys[[i]]),], na.rm=T)
    hetbaseline<-lm(hetresponse ~ 1)
    hetfit=lm(hetresponse ~ Yhats_het[[i]][rownames(Ys[[i]]),])
    
    het_test_stat.herit<--2*(logLik(hetbaseline))+2*logLik(hetfit)
    het_cv_pvals.herit[[i]]<-pchisq(het_test_stat.herit,1,lower.tail=F)
    het_cv_r2s.herit[[i]]<-summary(hetfit)$adj.r.squared
    
    
    m3<-lm(hom_expr_mat[rownames(Ys[[i]]),i] ~ Yhats_het[[i]][rownames(Ys[[i]]),])
    het_scales[[i]]=coef(m3)[2]
    # full model
    m4=lm(hom_expr_mat[rownames(Ys[[i]]),i] ~ Yhats_hom[[i]][rownames(Ys[[i]]),] + 
            Yhats_het[[i]][rownames(Ys[[i]]),])
    # test for signif of full model
    full_test_stat<--2*(logLik(m1))+2*logLik(m4)
    full_cv_pvals[[i]]<-pchisq(full_test_stat,2,lower.tail=F)
    full_cv_r2s[[i]]<-summary(m4)$adj.r.squared
    full_weights[[i]]<-coef(m4)[2:3]
    # test for signif of het | hom
    het_test_stat<--2*(logLik(m2))+2*logLik(m4)
    het_cv_pvals[[i]]<-pchisq(het_test_stat,1,lower.tail=F)
    het_cv_r2s[[i]]<-summary(m3)$adj.r.squared
    # test for signif of hom | het
    hom_test_stat<--2*(logLik(m3))+2*logLik(m4)
    hom_cv_pvals[[i]]=pchisq(hom_test_stat,1,lower.tail=F)
    hom_cv_r2s[[i]]=summary(m2)$adj.r.squared
  }else{
    # we didn't get a heterogeneous model
    het_cv_pvals[[i]]<-NA
    het_cv_r2s[[i]]<-NA
    het_cv_pvals.herit[[i]]<-NA
    het_cv_r2s.herit[[i]]<-NA
    # see if the homogeneoous explains some variability
    if(runcontent){
      hom_test_stat<--2*(logLik(m1))+2*logLik(m2)
      hom_cv_pvals[[i]]=pchisq(hom_test_stat,1,lower.tail=F)
      hom_cv_r2s[[i]]=summary(m2)$adj.r.squared
    }else{
      hom_cv_pvals=NA
      hom_cv_r2s=NA
    }
    full_cv_pvals[[i]]<-NA
    full_cv_r2s[[i]]<-NA
    full_weights[[i]]<-c(NA, NA)
  }
  # tissue by tissue approach
  if(!all(is.na(Yhats_tiss[[i]]))){
    t1<-lm(hom_expr_mat[rownames(Ys[[i]]),i] ~ Yhats_tiss[[i]][rownames(Ys[[i]]),])
    tiss_test_stat<--2*(logLik(m1))+2*logLik(t1)
    tiss_cv_pvals[[i]]<-pchisq(tiss_test_stat,1,lower.tail=F)
    tiss_cv_r2s[[i]]<-summary(t1)$adj.r.squared
  }
}
pvaldf=cbind(het_cv_pvals, het_cv_pvals.herit, hom_cv_pvals, hom_cv_pvals.herit, full_cv_pvals, tiss_cv_pvals)
rownames(pvaldf)=names(Ys)
r2df=cbind(het_cv_r2s, het_cv_r2s.herit, hom_cv_r2s, hom_cv_r2s.herit, full_cv_r2s, tiss_cv_r2s)
rownames(r2df)=names(Ys)
save(pvaldf, r2df, file = paste0(out_dir,gene_name,"_crossval_performance"))



if(verbose){
  message("Finished cross-validation, building TWAS weights.")
}

if(!dir.exists(paste0(twas_dir, "significant"))){
  system(paste0("mkdir ", twas_dir,"significant"))
}
# if(!dir.exists(paste0(twas_dir, "notsignificant"))){
#   system(paste0("mkdir ", twas_dir,"notsignificant"))
# }

# Now fit the model on all of the data for TWAS
## Fit the homogeneous component

if(runcontent){
  hom_fit<-big_spLinReg(X = explanatory, ind.train = match(rownames(hom_expr_mat), rownames(X)),
                        y.train=rowMeans(x=hom_expr_mat, na.rm = T), alphas=c(content_alpha),K=10,warn=F)
  hom_beta_vals<-unlist(summary(hom_fit)$beta[
    which.min(summary(hom_fit)$validation_loss)])[1:ncol(explanatory)]
  if(any(is.na(hom_beta_vals))){
    hom_beta_vals[which(is.na(hom_beta_vals))]<-0
  }
  if(length(hom_beta_vals)<dim(explanatory)[2]){
    new_betas<-rep(0, dim(explanatory)[2])
    new_betas[attr(hom_fit, "ind.col")]<-hom_beta_vals
    hom_beta_vals<-new_betas
  }
  ## Fit the heterogeneous components
  het_tiss_betas<-vector("list", q)
}

tiss_betas<-vector("list", q)
for(j in 1:q){
  if(runcontent){
    ### Find out which individuals were not present in other tissues
    nan_names<-names(rowMeans(hom_expr_mat[rownames(Ys[[j]]),-j],
                              na.rm = T)[which(is.nan(rowMeans(hom_expr_mat[rownames(Ys[[j]]),-j], na.rm = T)))])
    ### It may happen that there are individuals in this training set who
    ### do not appear in any other tissues, or they only appear in the test
    ### set of other tissues. If that's the case, remove them
    if(length(nan_names) == 0){
      subset_inds<-intersect(rownames(Ys[[j]]), rownames(hom_expr_mat))
      cur_hom_expr_mat<-hom_expr_mat[subset_inds,,drop=F]
    }else{
      subset_inds<-rownames(Ys[[j]])[!(rownames(Ys[[j]]) %in% nan_names)]
      subset_inds<-intersect(subset_inds, rownames(hom_expr_mat))
      cur_hom_expr_mat<-hom_expr_mat[subset_inds,,drop=F]
    }
    ### First, keep only the train individuals that are in this tissue+other tissues
    
    if(nrow(cur_hom_expr_mat) < 15){
      het_tiss_betas[[j]]<-rep(NA, m)
      ## Fit the tissue by tissue approach as well 
      tiss_response=hom_expr_mat[rownames(Ys[[j]]),j]
      
      tiss_fit<-big_spLinReg(X = explanatory, match(rownames(Ys[[j]]), rownames(X)),
                             y.train = tiss_response,K=10, alphas = c(tissue_alpha),warn=F)
      tiss_beta_vals<-unlist(summary(tiss_fit)$beta[
        which.min(summary(tiss_fit)$validation_loss)])[1:m]
      idx_remove<-c()
      if(any(is.na(tiss_beta_vals))){
        idx_remove<-which(is.na(tiss_beta_vals))
        tiss_beta_vals[idx_remove]<-0
      }
      if(length(tiss_beta_vals)<dim(explanatory)[2]){
        new_betas<-rep(0, dim(explanatory)[2])
        new_betas[attr(tiss_fit, "ind.col")]<-tiss_beta_vals
        tiss_beta_vals<-new_betas
      }
      tiss_betas[[j]]<-tiss_beta_vals
      next
    }
    
    ### tissue_j expression - mean(all tissues except j expression)
    het_response<-cur_hom_expr_mat[, j]-rowMeans(x=cur_hom_expr_mat,  na.rm = T)
    # het_explanatory<-X[rownames(cur_hom_expr_mat),]
    # het_explanatory<-as_FBM(het_explanatory, 
    #                         backingfile = paste0(out_dir,
    #                                              "content_tmp_",cur_fold))
    if(length(het_response) < 15){
      het_tiss_betas[[j]]<-rep(NA, m)
      ## Fit the tissue by tissue approach as well 
      tiss_response=hom_expr_mat[rownames(Ys[[j]]),j]
      
      tiss_fit<-big_spLinReg(X = explanatory, match(rownames(Ys[[j]]), rownames(X)),
                             y.train = tiss_response,K=10, alphas = c(tissue_alpha),warn=F)
      tiss_beta_vals<-unlist(summary(tiss_fit)$beta[
        which.min(summary(tiss_fit)$validation_loss)])[1:m]
      idx_remove<-c()
      if(any(is.na(tiss_beta_vals))){
        idx_remove<-which(is.na(tiss_beta_vals))
        tiss_beta_vals[idx_remove]<-0
      }
      if(length(tiss_beta_vals)<dim(explanatory)[2]){
        new_betas<-rep(0, dim(explanatory)[2])
        new_betas[attr(tiss_fit, "ind.col")]<-tiss_beta_vals
        tiss_beta_vals<-new_betas
      }
      tiss_betas[[j]]<-tiss_beta_vals
      next
    }
    
    het_fit<-big_spLinReg(X = explanatory, ind.train = match(rownames(cur_hom_expr_mat), rownames(X)),
                          y.train = het_response,K=10,alphas = c(content_alpha),warn=F)
    
    het_beta_vals<-unlist(summary(het_fit)$beta[
      which.min(summary(het_fit)$validation_loss)])[1:m]
    idx_remove<-c()
    if(any(is.na(het_beta_vals))){
      idx_remove<-which(is.na(het_beta_vals))
      het_beta_vals[idx_remove]<-0
    }
    if(length(het_beta_vals)<dim(explanatory)[2]){
      new_betas<-rep(0, dim(explanatory)[2])
      new_betas[attr(het_fit, "ind.col")]<-het_beta_vals
      het_beta_vals<-new_betas
    }
    het_tiss_betas[[j]]<-het_beta_vals
    
  }
  
  
  ## Fit the tissue by tissue approach as well 
  tiss_response=hom_expr_mat[rownames(Ys[[j]]),j]
  # tiss_explanatory=as_FBM(X[rownames(Ys[[j]]),],
  #                         backingfile = paste0(out_dir,
  #                                              "content_tmp_",cur_fold))
  tiss_fit<-big_spLinReg(X = explanatory, match(rownames(Ys[[j]]), rownames(X)),
                         y.train = tiss_response,K=10, alphas = c(tissue_alpha),warn=F)
  tiss_beta_vals<-unlist(summary(tiss_fit)$beta[
    which.min(summary(tiss_fit)$validation_loss)])[1:m]
  idx_remove<-c()
  if(any(is.na(tiss_beta_vals))){
    idx_remove<-which(is.na(tiss_beta_vals))
    tiss_beta_vals[idx_remove]<-0
  }
  if(length(tiss_beta_vals)<dim(explanatory)[2]){
    new_betas<-rep(0, dim(explanatory)[2])
    new_betas[attr(tiss_fit, "ind.col")]<-tiss_beta_vals
    tiss_beta_vals<-new_betas
  }
  tiss_betas[[j]]<-tiss_beta_vals
}

# Build CONTENT hom weights
if(runcontent){
  hom_weights<-hom_beta_vals
  snps=storesnps
  if(any(nchar(snps[,5])>1) | any(nchar(snps[,6])>1)){
    multichar<-unique(c(which(nchar(snps[,5])>1),
                        which(nchar(snps[,6])>1)))
    snps<-snps[-multichar,]
    hom_weights<-hom_weights[-multichar]
  }
  not_use<-which(is.na(snps[,2]))
  if(length(not_use)>0){
    hom_weights<-hom_weights[-not_use]
    snps<-snps[-not_use,]
  }
  wgt.matrix<-matrix(hom_weights, nrow=length(hom_weights), ncol=1)
  colnames(wgt.matrix)<-"CONTENT"
  rownames(wgt.matrix)=snps[,2]
  cv.performance=matrix(c(hom_cv_r2s.herit,hom_cv_pvals.herit), nrow=2, ncol=1, 
                        dimnames = list(c("rsq", "pval"), c("CONTENT")))
  if(any(unlist(hom_cv_pvals.herit) <= signif)){
    save(wgt.matrix, snps, cv.performance, file = 
           paste0(twas_dir, "significant/", gene_name, ".CONTENT.hom"))
  }else{
    # save(wgt.matrix, snps, cv.performance, file = 
    #       paste0(twas_dir, "notsignificant/", gene_name, ".CONTENT.hom"))
  }
  # Build CONTENT het weights and tissue by tissue weights
  hom.wgt.matrix=wgt.matrix
}


for(i in 1:q){
  if(runcontent){
    het_weights<-het_tiss_betas[[i]] 
  }
  tiss_weights<-tiss_betas[[i]]
  snps=storesnps
  if(any(nchar(snps[,5])>1) | any(nchar(snps[,6]) >1)){
    multichar<-unique(c(which(nchar(snps[,5])>1),
                        which(nchar(snps[,6])>1)))
    snps<-snps[-multichar,]
    if(runcontent){
      het_weights<-het_weights[-multichar] 
    }
    tiss_weights<-tiss_weights[-multichar]
  }
  not_use<-which(is.na(snps[,2]))
  if(length(not_use)>0){
    if(runcontent){
      het_weights<-het_weights[-not_use] 
    }
    tiss_weights<-tiss_weights[-not_use]
    snps<-snps[-not_use,]
  }
  # First, finish tissue by tissue weights
  wgt.matrix<-matrix(tiss_weights, nrow=length(tiss_weights), ncol=1)
  colnames(wgt.matrix)<-"tiss_spec"
  rownames(wgt.matrix)=snps[,2]
  cv.performance=matrix(c(tiss_cv_r2s[[i]],tiss_cv_pvals[[i]]), nrow=2, ncol=1, 
                        dimnames = list(c("rsq", "pval"), c("tiss_spec")))
  if(tiss_cv_pvals[[i]] <= signif){
    save(wgt.matrix, snps, cv.performance, file = 
           paste0(twas_dir, "significant/", gene_name, ".Tiss_spec.",names(Ys)[i]))
  }else{
    # save(wgt.matrix, snps, cv.performance, file = 
    #       paste0(twas_dir, "notsignificant/", gene_name, ".Tiss_spec.",names(Ys)[i]))
  }
  
  if(runcontent){
    wgt.matrix<-matrix(het_weights*het_scales[[i]], nrow=length(het_weights), ncol=1)
    colnames(wgt.matrix)<-"CONTENT"
    rownames(wgt.matrix)=snps[,2]
    cv.performance=matrix(c(het_cv_r2s.herit[[i]],het_cv_pvals.herit[[i]]), nrow=2, ncol=1, 
                          dimnames = list(c("rsq", "pval"), c("CONTENT")))
    if(!is.na(het_cv_pvals.herit[[i]])){
      if(het_cv_pvals.herit[[i]] <= signif){
        save(wgt.matrix, snps, cv.performance, file = 
               paste0(twas_dir, "significant/", gene_name, ".CONTENT_het.",names(Ys)[i]))
      }else{
        # save(wgt.matrix, snps, cv.performance, file = 
        #       paste0(twas_dir, "notsignificant/", gene_name, ".CONTENT_het.",names(Ys)[i]))
      }
      wgt.matrix<-matrix(het_weights, nrow=length(het_weights), ncol=1)
      wgt.matrix=wgt.matrix*full_weights[[i]][2] + hom.wgt.matrix*full_weights[[i]][1]
      colnames(wgt.matrix)<-"CONTENT"
      rownames(wgt.matrix)=snps[,2]
      cv.performance=matrix(c(full_cv_r2s[[i]],full_cv_pvals[[i]]), nrow=2, ncol=1, dimnames = list(c("rsq", "pval"), c("CONTENT")))
      if(full_cv_pvals[[i]] <= signif){
        save(wgt.matrix, snps, cv.performance, file = paste0(twas_dir, "significant/", gene_name, ".CONTENT_full.",names(Ys)[i]))
      }else{
        # save(wgt.matrix, snps, cv.performance, file = paste0(twas_dir, "notsignificant/", gene_name, ".CONTENT_full.",names(Ys)[i]))
      }
    } 
  }
  
}
system(paste0("rm ", paste0(out_dir,gene_name, "_content_tmp.bk")))
message("Done!")
if(runcontent){
  message("CONTENT discovered ", 
          sum(unlist(full_cv_pvals) <= signif, na.rm=T), " gene-tissue pair(s)")
  message(" and ",
          sum(unlist(het_cv_pvals) <= signif, na.rm=T), " heterogeneous component(s).") 
}
message("Tissue by tissue discovered ", 
        sum(unlist(tiss_cv_pvals) <= signif,na.rm=T), " gene-tissue pair(s)")
