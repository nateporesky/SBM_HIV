library(seqinr)
library(msa)
library(Biostrings)
library(stringdist)
library(stats)
library(ecodist)
#setwd("/gpfs/fs2/scratch/nporesky/SBM")
setwd("/Users/nateporesky/Library/CloudStorage/Box-Box/AndersonLabShared/Nathan_Poresky/SBM_Code")
dir()
print("this is for the stripped alignment part 1")
fastafile<- read.fasta(file = "HIV-1_env_gap_editedFull.fasta", seqtype = "AA", as.string = TRUE, set.attributes = FALSE)

mySeqs<- read.alignment(file = "HIV-1_env_gap_editedFull.fasta", format = "fasta")

print(table(sapply(mySeqs$seq, nchar)))

table(sapply(mySeqs$seq, nchar))

numChar<-sapply(mySeqs$seq, nchar)
pullMe1<-grep(829, numChar)#value changed!!
pullMe2<-grep(830, numChar)#value changed!!
pullMe3<-grep(831, numChar)#value changed!!
pullMe4<-grep(832, numChar)#value changed!!
pullMe5<-grep(833, numChar)#value changed!!
pullMe6<-grep(834, numChar)#value changed!!
pullMe7<-grep(835, numChar)#value changed!!
pullMe8<-grep(836, numChar)#value changed!!
pullMe9<-grep(837, numChar)#value changed!!
pullMe10<-grep(838, numChar)#value changed!!
pullMe11<-grep(839, numChar)#value changed!!
pullMe12<-grep(840, numChar)#value changed!!
pullMe13<-grep(841, numChar)#value changed!!
pullMe14<-grep(842, numChar)#value changed!!
pullMe15<-grep(843, numChar)#value changed!!
pullMe16<-grep(844, numChar)#value changed!!
pullMe17<-grep(845, numChar)#value changed!!
pullMe18<-grep(846, numChar)#value changed!!
pullMe19<-grep(847, numChar)#value changed!!
pullMe20<-grep(848, numChar)#value changed!!
pullMe21<-grep(849, numChar)#value changed!!
pullMe22<-grep(850, numChar)#value changed!!
pullMe23<-grep(851, numChar)#value changed!!
pullMe24<-grep(852, numChar)#value changed!!
pullMe25<-grep(853, numChar)#value changed!!
pullMe26<-grep(854, numChar)#value changed!!
pullMe27<-grep(855, numChar)#value changed!!
pullMe28<-grep(856, numChar)#value changed!!
pullMe29<-grep(857, numChar)#value changed!!
pullMe30<-grep(858, numChar)#value changed!!
pullMe31<-grep(859, numChar)#value changed!!
pullMe32<-grep(860, numChar)#value changed!!
pullMe33<-grep(861, numChar)#value changed!!
pullMe34<-grep(862, numChar)#value changed!!
pullMe35<-grep(863, numChar)#value changed!!
pullMe36<-grep(864, numChar)#value changed!!
pullMe37<-grep(865, numChar)#value changed!!
pullMe38<-grep(866, numChar)#value changed!!
pullMe39<-grep(867, numChar)#value changed!!
pullMe40<-grep(868, numChar)#value changed!!
pullMe41<-grep(869, numChar)#value changed!!
pullMe42<-grep(870, numChar)#value changed!!
pullMe43<-grep(871, numChar)#value changed!!
pullMe44<-grep(872, numChar)#value changed!!
pullMe45<-grep(873, numChar)#value changed!!
pullMe46<-grep(874, numChar)#value changed!!
pullMe47<-grep(875, numChar)#value changed!!
pullMe48<-grep(876, numChar)#value changed!!
pullMe49<-grep(877, numChar)#value changed!!
pullMe50<-grep(878, numChar)#value changed!!
pullMe51<-grep(879, numChar)#value changed!!
pullMe52<-grep(880, numChar)#value changed!!
pullMe53<-grep(881, numChar)#value changed!!
pullMe54<-grep(882, numChar)#value changed!!
pullMe55<-grep(883, numChar)#value changed!!
pullMe56<-grep(884, numChar)#value changed!!
pullMe57<-grep(885, numChar)#value changed!!
pullMe58<-grep(886, numChar)#value changed!!
pullMe59<-grep(887, numChar)#value changed!!
pullMe60<-grep(888, numChar)#value changed!!
pullMe61<-grep(889, numChar)#value changed!!
pullMeAll<-c(pullMe1, pullMe2, pullMe3, pullMe4, pullMe5, pullMe6, pullMe7, pullMe8, pullMe9, pullMe10, pullMe11, pullMe12, pullMe13, pullMe14, pullMe15, pullMe16, pullMe17, pullMe18, pullMe19, pullMe20, pullMe21, pullMe22, pullMe23, pullMe24, pullMe25, pullMe26, pullMe27, pullMe28, pullMe29, pullMe30, pullMe31, pullMe32, pullMe33, pullMe34, pullMe35, pullMe36, pullMe37, pullMe38, pullMe39, pullMe40, pullMe41, pullMe42, pullMe43, pullMe44, pullMe45, pullMe46, pullMe47, pullMe48, pullMe49, pullMe50, pullMe51, pullMe52, pullMe53, pullMe54, pullMe55, pullMe56, pullMe57, pullMe58, pullMe59, pullMe60, pullMe61)#Change these
pullNames<-names(fastafile)[pullMeAll]
QC1file<-fastafile[c(which(names(fastafile) %in% pullNames))]
write.fasta(as.list(QC1file), names=mySeqs$nam, file="QC1file.fasta")

namList<-list()
for (i in 1:length(QC1file)) {
  myName<-names(QC1file[i])
  aSeq<-QC1file[[i]]
  hasX<-which(strsplit(aSeq, "")[[1]]=="X")
  if(length(hasX)>0){
    namList<-c(namList,myName)}
}
QC2file<-QC1file[c(which(!(names(QC1file) %in% namList)))]
write.fasta(QC2file, names(QC2file), "QC2file.fasta", open = "w", nbchar = 60)

myClustalOAlignment <- readAAStringSet("QC2file.fasta", format = "fasta")
musAlign<-msa(myClustalOAlignment, "ClustalOmega")
writeXStringSet(unmasked(musAlign), file="QC2file_edit_strip_align.fasta")

timestamp()
save.image(file = "/gpfs/fs2/scratch/nporesky/SBM/FullList_edit_strip_msa.RData")

myAlign<-read.alignment(file = "QC2file_edit_strip_align.fasta", format = "fasta")
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
#save.image(file = "/gpfs/fs2/scratch/nporesky/SBM/FullList_edit_strip_align.RData")


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

write.csv(myMat, file="myEditFullDM.csv")
save.image(file = "/gpfs/fs2/scratch/nporesky/SBM/FullList_edit_strip_DM.RData")

timestamp()

myMDS<-cmdscale(myMat, k=2, eig=T)

myMDS$GOF

timestamp()

myDF <- as.data.frame(myMDS["points"])

timestamp()

write.csv(myDF, file = "myMDS.csv")

timestamp()

save.image(file = "/gpfs/fs2/scratch/nporesky/SBM/edit_FullList.RData")

