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

myAlign<-read.alignment(file = "fastafileSub1.fasta", format = "fasta")
table(sapply(myAlign$seq, nchar))

viruAln<-myAlign
for (i in 1:length(viruAln$nam)) {
  name<-viruAln$nam[[i]]
  seq<-viruAln$seq[[i]]
  assign(name, seq)
  #A epitope 	
  a<-substr(seq, 649, 678)#changed from myAlign table
  A<-paste0(a)
  assign(name, A)
  if(i == 1){
    namList<-name}
  else {namList<-c(namList,name)}
  
  
}

timestamp() 


myMat<-matrix(data=NA, nrow=length(namList), ncol=length(namList))
colnames(myMat)<-namList
rownames(myMat)<-namList
for (i in 1:length(namList)){
  name<-paste0(namList[i],"_vec")
  v1<-get(namList[i])
  for (j in 1:length(namList)){
    v2<-get(namList[j])
    dis<-stringdist(v1, v2, method="h")
    if(j==1){vec<-dis} else{vec<-c(vec, dis)}
    assign(name,vec)}
  if(j==length(namList)){myMat[, i]<-vec}
}

timestamp()

write.csv(myMat, file="myEditFullDM_P1.csv")

save.image(file = "/gpfs/fs2/scratch/nporesky/SBM/FullList_edit_strip_DM_P1.RData")

timestamp()

myMDS<-cmdscale(myMat, k=2, eig=T)

myMDS$GOF

myDF <- as.data.frame(myMDS["points"])

timestamp()

write.csv(myDF, file = "myMDSP1.csv")

