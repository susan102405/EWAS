### Functions for EWAS Pipeline-2

suppressWarnings(rm(rbint))
suppressWarnings(rm(export_results))
suppressWarnings(rm(publishFormat))
suppressWarnings(rm(splitAutosomal))
suppressWarnings(rm(sigResults))
suppressWarnings(rm(statsummary))

suppressWarnings(rm(f.RLM.par))
suppressWarnings(rm(f.LM.par))
suppressWarnings(rm(f.LM_RES.par))
suppressWarnings(rm(f.LM_RES_INT.par))
suppressWarnings(rm(f.LM_CAT.par))
suppressWarnings(rm(f.LOGISTIC.par))
suppressWarnings(rm(f.GEE_lm.par))
suppressWarnings(rm(f.GEE_logistic.par))
suppressWarnings(rm(f.LMM.par))


## Function to perform rank-based inverse normal transform (INT) - Blom transform
rbint <- function(u) 
{
    n <- length(u)
    r <- rank(u, ties.method = "average")
    out <- stats::qnorm((r - 0.375)/(n - 2 * k + 1))
    return(out)
}

## Function to export results
export_results <- function(modresults, NAMES_LIST, result_folder, rounddigit = rounddigit){
  modresults <- data.frame(modresults)
  if(NAMES_LIST$modelname == "lmCat") 
  {
    suppressPackageStartupMessages(library(stringi))
    
    colnames(modresults) = c("Estimates","StdErr", "Stat","Pvalue", 
                             "TypeII_SumSq", "TypeII_F", "TypeII_Pvalue",
                             "LSMEAN", "LSMEAN_SE","Sample_Size",
                             "beta_Min","beta_1stQuartile","beta_Median","beta_Mean","beta_3rdQuartile","beta_Max","beta_IQR","beta_SD",
                             "M_Min","M_1stQuartile","M_Median","M_Mean","M_3rdQuartile","M_Max","M_IQR","M_SD")
    rowNames <- rownames(modresults)
    rowNames <- gsub(".VAR", ".", rowNames)
    rowNames_split <- lapply(stri_split_fixed(stri_reverse(rowNames), ".", n = 2), stri_reverse)
    cpgNames <- sapply(rowNames_split, function(x) x[2])
    groupNames <- sapply(rowNames_split, function(x) x[1])
    modresults <- data.frame(CpG = cpgNames, Group = groupNames, modresults)
    rownames(modresults) <- NULL
    badTest <- which(rowSums(is.na(modresults)) == 26)
  } else {
    colnames(modresults) = c("Estimates","StdErr", "Stat","Pvalue",
                             "Sample_Size","beta_Min","beta_1stQuartile", "beta_Median","beta_Mean","beta_3rdQuartile","beta_Max","beta_IQR","beta_SD",
                             "M_Min","M_1stQuartile", "M_Median","M_Mean","M_3rdQuartile","M_Max","M_IQR","M_SD")
    modresults <- data.frame(CpG = rownames(modresults), modresults)
    badTest <- which(rowSums(is.na(modresults)) == 21)
  }
  
  if(length(badTest) > 0) 
  {
    message(paste0("The following CpGs were failed in the test and thus removed: ", paste0(unique(modresults$CpG[badTest]), collapse = ", ")))
    modresults <- modresults[-badTest, ]
  } else {
    message("All tests were successful without error!")
  }
  modresults$Sample_Size = as.integer(modresults$Sample_Size)
  modresults <- publishFormat(modresults, rounddigit = rounddigit)
  
  saveRDS(modresults, file = file.path(result_folder, 
                                    paste0(NAMES_LIST$cohortname, "_", NAMES_LIST$Year, "_", NAMES_LIST$VAR,"_",NAMES_LIST$modelname,"_",
                                           NAMES_LIST$datatype,"_",NAMES_LIST$cells,"_", NAMES_LIST$nPC,"PC_", NAMES_LIST$tag,"_",Sys.Date(),".RDS"))) 
  message("EWAS results exported!")
  return(modresults)
}

round_pad <- function(x, digits=0)
{
 format(round(x, digits), nsmall=digits)
}
                         
publishFormat<-function(res, rounddigit = 3){
  est = as.numeric(res[,"Estimates"])
  se = as.numeric(res[,"StdErr"])
  p = as.numeric(res[,"Pvalue"])
  res$lower <- est-1.96*se
  res$upper <- est+1.96*se
  res$beta <- round_pad(est, rounddigit)
  res$CI <- paste0("(",round_pad(res$lower, rounddigit),", ",round_pad(res$upper, rounddigit),")")
  res$p <- round_pad(p, rounddigit+5)
  return(res)
}

splitAutosomal <- function(res, annot){
  cpg_auto <- as.character(annot$Name[!annot$chr %in% c("chrX", "chrY")])
  cpg_X <- as.character(annot$Name[annot$chr %in% c("chrX")])
  cpg_Y <- as.character(annot$Name[annot$chr %in% c("chrY")])
  
  if(length(cpg_auto) == 0) {
    message("No autosomal CpG found!")
    results_auto <- NA
  } else {
    results_auto <- res[which(res$CpG %in% cpg_auto),]
    message(paste0(length(unique(results_auto$CpG)), " autosomal CpGs."))
  }
  
  if(length(cpg_X) == 0) {
    message("No ChrX CpG found!")
    results_X <- NA
  } else {
    results_X <- res[which(res$CpG %in% cpg_X),]
    message(paste0(length(unique(results_X$CpG)), " X-chromosome CpGs."))
  }
  
  if(length(cpg_Y) == 0) {
    message("No ChrY CpG found!")
    results_Y <- NA
  } else {
    results_Y <- res[which(res$CpG %in% cpg_Y),]
    message(paste0(length(unique(results_Y$CpG)), " Y-chromosome CpGs."))
  }
  
  return(list(auto = results_auto, X = results_X, Y = results_Y))
}

sigResults <- function(results, annotcord, NAMES_LIST, psigcut = psigcut, rounddigit = rounddigit, qval = TRUE){
  if(is.null(nrow(results))) 
  {
    message("No result input found!")
    return(NULL)
  }
    
  if(NAMES_LIST$modelname == "lmCat")
  {
    results$p.FDR<-p.adjust(results$TypeII_Pvalue,"fdr")
    if(qval) results$qvalue<-qvalue(results$TypeII_Pvalue)$qvalues else results$qvalue <- NA 
    sigCpG <- na.omit(results[, c("CpG", "TypeII_Pvalue")])
    sigCpG <- subset(sigCpG, TypeII_Pvalue<psigcut)
    sigCpG <- as.character(sigCpG$CpG[order(sigCpG$TypeII_Pvalue)])
    mind <- match(results$CpG, sigCpG)
    results <- results[intersect(order(mind), which(!is.na(mind))), ]
    # Add annotation
    results = cbind(results,annotcord[match(results$CpG,annotcord$Name),])
  } else {
    results$p.FDR<-p.adjust(results$Pvalue,"fdr")
    if(qval) results$qvalue<-qvalue(results$Pvalue)$qvalues else results$qvalue <- NA 
    results<-results[results$Pvalue<psigcut,]
    results<-results[order(results$Pvalue),]
    # Add annotation
    results = cbind(results,annotcord[match(results$CpG,annotcord$Name),])
  }
  
  sigFileName <- paste0(NAMES_LIST$cohortname, "_", NAMES_LIST$Year, "_", NAMES_LIST$VAR,"_", NAMES_LIST$modelname,"_",
                        NAMES_LIST$datatype,"_", NAMES_LIST$cells,"_", NAMES_LIST$nPC,"PC_", NAMES_LIST$tag,"_",Sys.Date(),".csv")
  
  message(paste0("Writing file: ", sigFileName))
  write.csv(results, file.path(result_folder, sigFileName), row.names = FALSE)
  message("Signficant results exported!")
}

# Add summary of statistics of the tested CpG sites (will add 17 columns)
statsummary <- function(bigdata, type){
  samplesize <- nrow(bigdata)
  if(type == "Mval")
  {
    Mval <- bigdata$methy
    betaVal <- 2^Mval/(2^Mval + 1)
    res = c(samplesize, min(betaVal),quantile(betaVal,0.25, na.rm = TRUE),median(betaVal),mean(betaVal),quantile(betaVal,0.75, na.rm = TRUE),max(betaVal),IQR(betaVal),sd(betaVal),
          min(Mval),quantile(Mval,0.25, na.rm = TRUE),median(Mval),mean(Mval),quantile(Mval,0.75, na.rm = TRUE),max(Mval),IQR(Mval),sd(Mval))
    return(res)
  }

  if(type == "beta")
  {
    betaVal <- bigdata$methy
    Mval <- log2(betaVal/(1-betaVal))
    res = c(samplesize, min(betaVal),quantile(betaVal,0.25, na.rm = TRUE),median(betaVal),mean(betaVal),quantile(betaVal,0.75, na.rm = TRUE),max(betaVal),IQR(betaVal),sd(betaVal),
          min(Mval),quantile(Mval,0.25, na.rm = TRUE),median(Mval),mean(Mval),quantile(Mval,0.75, na.rm = TRUE),max(Mval),IQR(Mval),sd(Mval))
    return(res)
  }
    
  if(type == "Mval_residual")
  {
    resid <- bigdata$methy
    res = c(samplesize, min(resid),quantile(resid,0.25, na.rm = TRUE),median(resid),mean(resid),quantile(resid,0.75, na.rm = TRUE),max(resid),IQR(resid),sd(resid))
    return(res)
  }
    
  if(type == "beta_residual")
  {
    resid <- bigdata$methy
    res = c(samplesize, min(resid),quantile(resid,0.25, na.rm = TRUE),median(resid),mean(resid),quantile(resid,0.75, na.rm = TRUE),max(resid),IQR(resid),sd(resid))
    return(res)
  }      
}

### Modeling functions:

## Debug
# methcol = setNames(seq_len(ncol(tdatRUN)), dimnames(tdatRUN)[[2]])[1]

## RLM
f.RLM.par <- function(methcol, VAR, COV, model_statement, datatype, tdatRUN) { 
  bigdata <- data.frame(na.omit(cbind(VAR = eval(parse(text = paste0("df$", VAR))),methy = tdatRUN[, methcol], COV)))
  mod <- try(rlm(model_statement, bigdata, maxit=200))
  # pull out a data.frame with results
  if("try-error" %in% class(mod)){
    b <- rep(NA, 21)
  } else {
    cf <- try(coeftest(mod, vcov=vcovHC(mod, type="HC0")))
    if(class(cf) == "try-error"){
      b <- rep(NA, 21)
    } else {b <- c(cf[2,], statsummary(bigdata, datatype))
    }
  }
  invisible(b)
}

## LM
f.LM.par <- function(methcol, VAR, COV, model_statement, datatype, tdatRUN) { 
  bigdata <- data.frame(na.omit(cbind(VAR = eval(parse(text = paste0("df$", VAR))),methy = tdatRUN[, methcol], COV)))
  mod <- try(lm(model_statement, bigdata))
  if("try-error" %in% class(mod)){
    b <- rep(NA, 21)
  } else {
    cf <- summary(mod)$coefficients
    b <- c(cf[2,], statsummary(bigdata, datatype))
  }
  invisible(b)
}

## LM_RES
f.LM_RES.par <- function(methcol, VAR, COV, model_statement, res_model_statement, datatype, tdatRUN) { 
  bigdata <- data.frame(na.omit(cbind(VAR = eval(parse(text = paste0("df$", VAR))), methy = tdatRUN[, methcol], COV)))
  mod_res <- try(lm(res_model_statement, bigdata))
  if("try-error" %in% class(mod_res)){
    b <- rep(NA, 21)
  } else {
    bigdata$methy <- residuals(mod_res)
    mod <- try(lm(model_statement, bigdata))
    if("try-error" %in% class(mod)){
      b <- rep(NA, 21)
    } else {
      cf <- summary(mod)$coefficients
      b <- c(cf[2,], nrow(bigdata), rep(NA, 16)) # statsummary is not applicable for residual
    }
  }
  invisible(b)
}
    
## LM_RES_INT
f.LM_RES_INT.par <- function(methcol, VAR, COV, model_statement, res_model_statement, datatype, tdatRUN) { 
  bigdata <- data.frame(na.omit(cbind(VAR = eval(parse(text = paste0("df$", VAR))), methy = tdatRUN[, methcol], COV)))
  mod_res <- try(lm(res_model_statement, bigdata))
  if("try-error" %in% class(mod_res)){
    b <- rep(NA, 21)
  } else {
    bigdata$methy <- rbint(residuals(mod_res))
    mod <- try(lm(model_statement, bigdata))
    if("try-error" %in% class(mod)){
      b <- rep(NA, 21)
    } else {
      cf <- summary(mod)$coefficients
      b <- c(cf[2,], nrow(bigdata), rep(NA, 16)) # statsummary is not applicable for residual
    }
  }
  invisible(b)
}
                         
## LM_CAT
f.LM_CAT.par <- function(methcol, VAR, nCat, COV, model_statement, datatype, tdatRUN) { 
  bigdata <- data.frame(na.omit(cbind(VAR = eval(parse(text = paste0("df$", VAR))), methy = tdatRUN[, methcol], COV)))
  mod <- try(lm(model_statement, bigdata))
  if("try-error" %in% class(mod)){
    b <- rep(NA, 26)
  } else {
    anova_typeII <- as.matrix(Anova(mod))[1,]
    anova_typeII <- rbind(anova_typeII, matrix(rep(rep(NA, 4), nCat - 1), nrow = nCat-1))
    lsm <- emmeans(mod, ~VAR, data = bigdata)
    cf <- summary(mod)$coefficients[seq_len(nCat-1) + 1,]
    cf <- rbind(NA, cf)
    cf <- cbind(cf, anova_typeII[,-2], as.data.frame(summary(lsm))[,c("emmean", "SE")])
    statsummary_res <- NULL
    for(l in lsm@levels$VAR)
    {
      statsummary_res <- rbind(statsummary_res, statsummary(subset(bigdata, VAR == l), datatype))
    }
    b <- cbind(cf, statsummary_res)
  }
  invisible(b)
}

## LOGISTIC
f.LOGISTIC.par <- function(methcol, VAR, COV, model_statement, datatype, tdatRUN) { 
  bigdata <- data.frame(na.omit(cbind(VAR = eval(parse(text = paste0("df$", VAR))), methy = tdatRUN[, methcol], COV)))
  mod <- try(glm(model_statement, bigdata, family = binomial))
  if("try-error" %in% class(mod)){
    b <- rep(NA, 21)
  } else {
    cf <- summary(mod)$coefficients
    b <- c(cf[2,], statsummary(bigdata, datatype))
  }
  invisible(b)
}

## GEE-linear
f.GEE_LM.par <- function(methcol, VAR, COV, ID, model_statement, datatype, tdatRUN) { 
  bigdata <- data.frame(na.omit(cbind(VAR = eval(parse(text = paste0("df$", VAR))), methy = tdatRUN[, methcol], COV, ID = ID)))
  mod <- try(geeglm(model_statement, id = ID, data = bigdata, family = gaussian, corstr="ar1"))
  if("try-error" %in% class(mod)){
    b <- rep(NA, 21)
  } else {
    cf <- summary(mod)$coefficients
    b <- c(cf[2,], statsummary(bigdata, datatype))
  }
  invisible(b)
}

## GEE-logistic
f.GEE_LOGISTIC.par <- function(methcol, VAR, COV, ID, model_statement, datatype, tdatRUN) { 
  bigdata <- data.frame(na.omit(cbind(VAR = eval(parse(text = paste0("df$", VAR))), methy = tdatRUN[, methcol], COV, ID = ID)))
  mod <- try(geeglm(model_statement, id = ID, data = bigdata, family = binomial, corstr="ar1"))
  if("try-error" %in% class(mod)){
    b <- rep(NA, 21)
  } else {
    cf <- summary(mod)$coefficients
    b <- c(cf[2,], statsummary(bigdata, datatype))
  }
  invisible(b)
}

## LMM-linear mixed model
f.LMM.par <- function(methcol, VAR, COV, ID, model_statement, datatype, tdatRUN) { 
  bigdata <- data.frame(na.omit(cbind(VAR = eval(parse(text = paste0("df$", VAR))), methy = tdatRUN[, methcol], COV, ID = ID)))  
  mod <- try(lmer(model_statement, data = bigdata))
  if("try-error" %in% class(mod)){
    b <- rep(NA, 21)
  } else {    
    cf <- summary(mod)$coefficients    
    b <- c(cf[2,], statsummary(bigdata, datatype))
  }
  invisible(b)
}

message("Function2.R loaded!")
