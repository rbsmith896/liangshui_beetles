#RB Smith
#MSc Thesis 2024: Liangshui Beetles
#Original Functions

library(dplyr)
library(rarestR)
library(ggplot2)
library(abdiv)
library(tidyverse)
library(viridis)
library(fossil)
library(vegan)
library(scales)
library(gridExtra)

#calculate Chao1 index for all plots in a site x species matrix
getChaos <- function(beets_mat){
  output <- matrix(0,nrow(beets_mat),1)
  row.names(output) <- row.names(beets_mat)
  for(i in 1:nrow(output)){
    output[i,1] <- chao1(beets_mat[i,])
  }
  return(output)
}

#Get TESa, b, and ab values for all plots in a dataframe
TESify <- function(beets_df){
  p <- dim(beets_df)[1]
  output <- data.frame(TESa.est = double(p),TESa.sd = double(p),TESa.model = character(p),
                       TESb.est = double(p),TESb.sd = double(p),TESb.model = character(p),
                       TESab.est = double(p),TESab.sd = double(p),
                       alpha = double(p), TESa.cov = double(p), TESb.cov = double(p), TESab.cov = double(p))
  
  row.names(output) <- row.names(beets_df)
  output$alpha <- rowSums(beets_df != 0)
  
  beets_mat <- as.matrix(beets_df)
  
  for(i in 1:p){
    tes_out <- tes(beets_mat[i,])
    
    output$TESa.est[i] <- tes_out$tbl[1,1]
    output$TESa.sd[i] <- tes_out$tbl[1,2]
    output$TESa.model[i] <- tes_out$tbl[1,3]
    
    output$TESb.est[i] <- tes_out$tbl[2,1]
    output$TESb.sd[i] <- tes_out$tbl[2,2]
    output$TESb.model[i] <- tes_out$tbl[2,3]
    
    output$TESab.est[i] <- tes_out$tbl[3,1]
    output$TESab.sd[i] <- tes_out$tbl[3,2]
    
    output$TESa.cov[i] <- output$alpha[i]/tes_out$tbl[1,1]
    output$TESb.cov[i] <- output$alpha[i]/tes_out$tbl[2,1]
    output$TESab.cov[i] <- output$alpha[i]/tes_out$tbl[3,1]
  }
  
  output$sampleSize <- as.vector(rowSums(beets_df))
  return(output)
}

#Get a pairwise matrix of TESS species shared estimates
TESSdist <- function(beets_df){
  output <- ess(beets_df, m=1, index = "CNESS")
  beets_mat <- as.matrix(beets_df)
  count <- 1
  print(nrow(beets))
  for(i in 1:(nrow(beets_df)-1)){
    for(j in (i+1):nrow(beets_df)){
      print("Combo:")
      print(row.names(beets_df)[i])
      print(row.names(beets_df)[j])
      res <- try(tess(rbind(beets_mat[i,],beets_mat[j,])))
      if(inherits(res, "try-error")){
        output[count] <- NA
        count <- count+1
      }
      else{
        tess_out <- tess(rbind(beets_mat[i,],beets_mat[j,]))
        output[count] <- as.double(tess_out$tbl[1])
        count <- count+1
      }
    }
  }
  return(output)
}

#take the TESS dist output and turn it into a literal matrix
makeTESSMatrix <- function(TESSDist, plots){
  tess_matrix <- matrix(0,length(plots),length(plots))
  row.names(tess_matrix) <- plots
  colnames(tess_matrix) <- plots
  count <- 1
  for(i in 1:(length(plots)-1)){
    for(j in (i+1):length(plots)){
      tess_matrix[j,i] <- TESSDist[count]
      count <- count+1
    }
  }
  return(tess_matrix)
}

#now count the number of NA results for each plot
countNAs <- function(tess_mat){
  output <- matrix(0,nrow(tess_mat),1)
  row.names(output) <- row.names(tess_mat)
  for(i in 1:nrow(output)){
    output[i,1] <- sum(is.na(tess_mat[i,])) + sum(is.na(tess_mat[,i]))
  }
  return(output)
}

#Get a TESS matrix without NAs by recursively removing problematic plots
#She's slow but she'll get there
getWorkingTESSdist <- function(beets_df){
  tessdistmatrix <- TESSdist(beets_df)
  tess_mat <- makeTESSMatrix(tessdistmatrix,rownames(beets_df))
  NAcounts <- countNAs(tess_mat)
  count <- 0
  rerun.beets <- beets_df
  rerun.plotnames <- rownames(beets_df)
  rerun.NAcounts <- NAcounts
  while(sum(rerun.NAcounts!=0)){
    print("=============================================================================")
    print("=============================================================================")
    print(paste("RERUN NUMBER: ",count,sep=""))
    print("=============================================================================")
    print("=============================================================================")
    rerun.beets <- rerun.beets[-which.max(rerun.NAcounts[,1]),]
    rerun.plotnames <- rerun.plotnames[-which.max(rerun.NAcounts[,1])]
    rerun.tessdistmatrix <- TESSdist(rerun.beets)
    rerun.tess_mat <- makeTESSMatrix(rerun.tessdistmatrix,rerun.plotnames)
    rerun.NAcounts <- countNAs(rerun.tess_mat)
    count <- count+1
    if(count>nrow(beets_df)){
      print("FUNCTION TIMED OUT")
      NAcounts <- c(0)
    }
  }
  return(rerun.tess_mat)
}
