install.packages("BiocManager")
library(BiocManager)
install(c("genbankr","Biostrings","ggtree","annotate")) 
install.packages(c("reshape2","rentrez"))
library(genbankr)
#find ncov in Gen bank using succession number
ncovID<-GBAccession("MN994468")
ncovID
#read gen bank to find ID
ncovk<-readGenBank(ncovID)
ncovk
str(ncovk)
class(ncovk)
otherFeatures(ncovk)
library(annotate)
#sequence blast
ncovkBLAST<-blastSequences(paste(ncovk@sequence),as = 'data.frame',hitListSize = 40, timeout = 60)
ncovkBLAST
#find hits in query
ncovHitsDF<-data.frame(ID=ncovkBLAST$Hit_accession,Seq=ncovkBLAST$Hsp_hseq,stringsAsFactors = F)
ncovHitsDF
ncovkBLAST$Hit_len
#read hit sequences in query
ncovHitSeqs<-read.GenBank(ncovkBLAST$Hit_accession[1:3])
ncovHitSeqs
attr(ncovHitSeqs,"species")
head(ncovHitSeqs)
head(ncovkBLAST)
# NJ phylogenetic tree for unknown (unk)
library(ggtree)
install.packages("ape")
library(ape)
as.character(ncovHitSeqs)
length(ncovHitSeqs$CP002746)
ncovHitsDF[ncovHitsDF$ID=="CP002746",]
#new object with separate columns and split DNA sequence
ncovHitsDNA<-sapply(ncovHitsDF$Seq,strsplit,split="")
ncovHitsDNA
#give sequence a name
names(ncovHitsDNA)<-paste(1:nrow(ncovHitsDF),ncovHitsDF$ID,sep="_")
#conver to DNAbin
ncovHitsDNA<-as.DNAbin(ncovHitsDNA)
#run muscle
ncovAlign<-muscle(ncovHitsDNA,quiet=F)
#checkalignment
checkAlignment(ncovAlign[1:20,1:100],what=1)
checkAlignment(ncovAlign,what=3)
KeepSeq<-SeqLen>1000
ncovSubset<-ncovAlign[KeepSeq,]
checkAlignment(ncovSubset,what=1)
ncovSubAlign<-muscle(ncovSubset,quiet=F)
#distance matrix
ncovDM<-dist.dna(ncovSubAlign,model="K80")
class(ncovDM)
length(ncovDM)
#generating tree
ncovTree<-nj(ncovDM)
str(ncovTree)
class(ncovTree)
library(ggtree)
ggtree(ncovTree)
