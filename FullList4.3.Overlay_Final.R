library(seqinr)
library(msa)
library(Biostrings)
library(stringdist)
library(stats)
library(ecodist)
#setwd("/gpfs/fs2/scratch/nporesky/SBM")
setwd("/Users/nateporesky/Library/CloudStorage/Box-Box/AndersonLabShared/Nathan_Poresky/SBM_Code")
dir()

#Randomizing Fasta Files
fastafile<- read.fasta(file = "QC2file_edit_strip_align.fasta", seqtype = "AA", as.string = TRUE, set.attributes = FALSE)
fastafile_shuffle <- fastafile[sample(length(fastafile))]
fastafileSub1<-fastafile_shuffle[11:35010] 
fastafileSub2<-fastafile_shuffle[35011:70010]
fastafileSub3<-fastafile_shuffle[70011:95469]
write.fasta(as.list(fastafileSub1), names=names(fastafileSub1), file="fastafileSub1.0.fasta")
write.fasta(as.list(fastafileSub2), names=names(fastafileSub2), file="fastafileSub2.0.fasta")
write.fasta(as.list(fastafileSub3), names=names(fastafileSub3), file="fastafileSub3.0.fasta")
shared_sequences <- fastafile_shuffle[1:10]
write.fasta(as.list(shared_sequences), names=names(shared_sequences), file="shared_sequences.fasta")


##Creating CMDScale
#Group1
myAlign<-read.alignment(file = "fastafileSub1.fasta", format = "fasta")
table(sapply(myAlign$seq, nchar))

viruAln<-myAlign
for (i in 1:length(viruAln$nam)) {
  name<-viruAln$nam[[i]]
  seq<-viruAln$seq[[i]]
  assign(name, seq)
  #A epitope 	
  a<-substr(seq, 1, 1460)#changed from myAlign table
  A<-paste0(a)
  assign(name, A)
  if(i == 1){
    namList<-name}
  else {namList<-c(namList,name)}
  
  
}


timestamp() 
save.image(file = "/gpfs/fs2/scratch/nporesky/SBM/FullList_edit_strip_align.RData")



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

timestamp()

myDF <- as.data.frame(myMDS["points"])

timestamp()

write.csv(myDF, file = "myMDSP1.csv")



#Group2
myAlign<-read.alignment(file = "fastafileSub2.fasta", format = "fasta")
table(sapply(myAlign$seq, nchar))

viruAln<-myAlign
for (i in 1:length(viruAln$nam)) {
  name<-viruAln$nam[[i]]
  seq<-viruAln$seq[[i]]
  assign(name, seq)
  #A epitope 	
  a<-substr(seq, 1, 1460)#changed from myAlign table
  A<-paste0(a)
  assign(name, A)
  if(i == 1){
    namList<-name}
  else {namList<-c(namList,name)}
  
  
}


timestamp() 
save.image(file = "/gpfs/fs2/scratch/nporesky/SBM/FullList_edit_strip_align.RData")



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

write.csv(myMat, file="myEditFullDM_P2.csv")

save.image(file = "/gpfs/fs2/scratch/nporesky/SBM/FullList_edit_strip_DM_P2.RData")

timestamp()

myMDS<-cmdscale(myMat, k=2, eig=T)

myMDS$GOF

timestamp()

myDF <- as.data.frame(myMDS["points"])

timestamp()

write.csv(myDF, file = "myMDSP2.csv")





#Group 3
myAlign<-read.alignment(file = "fastafileSub3.fasta", format = "fasta")
table(sapply(myAlign$seq, nchar))

viruAln<-myAlign
for (i in 1:length(viruAln$nam)) {
  name<-viruAln$nam[[i]]
  seq<-viruAln$seq[[i]]
  assign(name, seq)
  #A epitope 	
  a<-substr(seq, 1, 1460)#changed from myAlign table
  A<-paste0(a)
  assign(name, A)
  if(i == 1){
    namList<-name}
  else {namList<-c(namList,name)}
  
  
}


timestamp() 
save.image(file = "/gpfs/fs2/scratch/nporesky/SBM/FullList_edit_strip_align.RData")



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

write.csv(myMat, file="myEditFullDM_P3.csv")

save.image(file = "/gpfs/fs2/scratch/nporesky/SBM/FullList_edit_strip_DM_P3.RData")

timestamp()

myMDS<-cmdscale(myMat, k=2, eig=T)

myMDS$GOF

timestamp()

myDF <- as.data.frame(myMDS["points"])

timestamp()

write.csv(myDF, file = "myMDSP3.csv")


#Plotting


myMDSP1 <- read.csv("myMDSP1.csv")
plot(myMDSP1$points.1,myMDSP1$points.2)

myMDSP1.5 <- read.csv("myMDSP1.5.csv")
plot(myMDSP1$points.1,myMDSP1$points.2)

myMDSP2 <- read.csv("myMDSP2.csv")
plot(myMDSP2$points.1,myMDSP2$points.2)

myMDSP3 <- read.csv("myMDSP3.csv")
plot(myMDSP3$points.1,myMDSP3$points.2)

