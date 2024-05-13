library(seqinr)
library(msa)
library(Biostrings)
library(stringdist)
library(stats)
library(ecodist)
setwd("/gpfs/fs2/scratch/nporesky/SBM")
#setwd("/Users/nateporesky/Library/CloudStorage/Box-Box/AndersonLabShared/Nathan_Poresky/SBM_Code")

dir()
print("this is for the Group 1 sequences for Data Matrix")

load(file = "/gpfs/fs2/scratch/nporesky/SBM/FullList_edit_strip_DM_P1.RData")
timestamp()

myMDS<-cmdscale(myMat, k=2, eig=T)

myMDS$GOF

myDF <- as.data.frame(myMDS["points"])

timestamp()

write.csv(myDF, file = "myMDSP1.csv")

