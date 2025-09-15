#!/usr/bin/env Rscript  
library(zeallot)
library(multtest)
library(gplots)
library(LDheatmap)
library(genetics)
library(ape)
library(EMMREML)
library(compiler)
library("scatterplot3d")

c(baseDir, mode, traits, hmp, pcaDataCSV=NULL)%<-% commandArgs(TRUE)

source(paste(baseDir, "gapit_functions.txt_new", sep="/"))
source(paste(baseDir, "emma.txt_new", sep="/"))
#source(paste(baseDir, "FarmCPU_functions.txt", sep="/"))

myY <-read.table(traits, head = TRUE)
myG <-read.delim(hmp, head = FALSE)
mode1 <- c(mode)

if (mode == "FASTLMM") {
  if (!is.null(pcaDataCSV)) {
    for(i in 2:ncol(myY))
  {
    myCV <-read.csv(pcaDataCSV, head=TRUE, sep=",")
    myGAPIT_SUPER <-GAPIT(
    #Y=myY[,c(1,2)],
    Y=myY[,c(1,i)],
    G=myG,
    CV=myCV,
    sangwich.top="SUPER", #options are GLM,MLM,CMLM, FaST and SUPER 
    sangwich.bottom="SUPER", #options are GLM,MLM,CMLM, FaST and SUPER
    LD=0.1)
    }
  } else if (is.null(pcaDataCSV)) {
    for(i in 2:ncol(myY))
  {
    myGAPIT_SUPER <-GAPIT(
   #Y=myY[,c(1,2)],
    Y=myY[,c(1,i)],
    G=myG,
    sangwich.top="SUPER", #options are GLM,MLM,CMLM, FaST and SUPER 
    sangwich.bottom="SUPER", #options are GLM,MLM,CMLM, FaST and SUPER
    LD=0.1)
  }
  }
} else if (mode == "BICmodelSelection") {
  myGAPIT <- GAPIT(
  Y=myY,
  G=myG,
  PCA.total=10,
  Model.selection = TRUE
  )
} else if (mode != "FASTLMM") {
  mode2 <- unlist(strsplit(mode1,","))
  for(i in 2:ncol(myY))
  {
  myGAPIT <- GAPIT(
  Y=myY[,c(1,i)],
  G=myG,
  PCA.total=3,
  #Multiple_analysis = TRUE,
  model = c(mode2)#,
  #method.bin="optimum",
  #cutOff = 0.05,
  #bin.size=c(5e5,5e6,5e7),
  #bin.selection=seq(10,100,10)
  )
  }
}

#else if (mode == "cMLM") {
  #myGAPIT <- GAPIT(
  #Y=myY,
  #G=myG,
  #PCA.total=3
  #)
#} else if (mode == "GLM") {
#  myGAPIT <- GAPIT(
#  Y=myY,
 # G=myG,
 # PCA.total=3,
 # model="GLM"
 # )
#}

