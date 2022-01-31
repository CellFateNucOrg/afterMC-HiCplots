# create table of consecutive seqContacts

#library(rhdf5)
### Initiation variables, packages, etc
library(lattice)
require(stats)
library(tibble)
library(dplyr)
library(tidyr)
library(ggplot2)
library(gridExtra)
library(GenomicRanges)
library(BSgenome.Celegans.UCSC.ce11)
library(viridis)
library(ggpubr)

source("./js_fragPlottingFunctions.R")



options(scipen=1e8)
options(tibble.width = Inf)
#options(digits=9)
Chrnames<-c("chrI","chrII","chrIII","chrIV","chrV","chrX") # used to get rid of mtDNA
wbChr<-c("I","II","III","IV","V","X","MtDNA") # used to reorder sequence name levels in objects
importantStats<-list()
dir.create("./plots_js")
dir.create("./rds")
dir.create("./csv")


sampleName="20190501_HIC6_7_barcode08_pass"
path="./"

importantStats<-readRDS(paste0("./rds/importantStats_",sampleName,".rds"))


##############################################################################
### make table of sequential seqContacts in read (seqContacts)
#############################################################################


###############
##############
###############

### hop table

# hop table has several additional fields as follows:
# AlnChr2 - chromosome to which the second fragment aligns to
# hopType - within same chromosome (Intra) or between chromosomes (Inter)
# hopSize - minimal absolute distance between fragment ends
# sameStrand - are both fragments on the both Strand?
# pairOri - orientation of hop:
#       tandem:     -frag1-->  -frag2-->  or  <--frag2- <--frag1-
#       inverted:   -frag2--> -frag1-->   or  <--frag1- <--frag2-
#       convergent: -frag1--> <--frag2-   or  -frag2--> <--frag1-
#       divergent:  <--frag2- -frag1-->   or  <--frag1- -frag2-->
# hopOverlap - do the aligned fragments overlap?
# firstChr - chr,Strand,start and end of alignment of first fragment in pair
# secondChr - chr, Strand, start and end of alignment of second fragment in pair

FragAlin<-readRDS(paste0("./rds/multiMapped_",sampleName,".rds"))

#FragAlin<-read.csv(paste0(path, expName,"_realFrag.csv"),header=T,stringsAsFactors=F)

seqContacts<-FragAlin %>%
  distinct(ReadID, FrgID, .keep_all = TRUE) %>% # remove duplicated fragments
  arrange(ReadID,FrgID) %>% # sort by read id then fragment id
  mutate(RefStart1=RefStart,
         RefEnd1=RefEnd,
         RefStart2=lead(RefStart),
         RefEnd2=lead(RefEnd),
         RefStrand1=Strand,
         RefStrand2=lead(Strand),
         RefChr1=Chr,
         RefChr2=lead(Chr)) %>%
  filter(!is.na(RefStart2))

seqContacts$Chr<-NULL
seqContacts$RefStart<-NULL
seqContacts$RefEnd<-NULL
seqContacts$Strand<-NULL


seqContacts$pairOri<-pairOrientation(with(seqContacts,cbind(RefStart1,RefEnd1,
                                         RefStart2,RefEnd2,
                                         RefStrand1,RefStrand2,
                                         RefChr1,RefChr2)))

seqContacts$overlap<-isOverlapping(with(seqContacts,cbind(RefStart1,RefEnd1,
                                                      RefStart2,RefEnd2,
                                                      RefStrand1,RefStrand2,
                                                      RefChr1,RefChr2)))
seqContacts<-seqContacts %>%
  mutate(hopType=ifelse(RefChr1==RefChr2,"cis","trans")) %>%
  mutate(hopSize=minDistance(RefStart1,RefEnd1,RefStart2,RefEnd2)) %>%
  mutate(sameStrand=ifelse(RefStrand1==RefStrand2,TRUE,FALSE)) %>%
  mutate(firstChr=paste0(wbChr[RefChr1],ifelse(RefStrand1=="True","+","-"),
                         ":",RefStart1,"-",RefEnd1)) %>%
  mutate(secondChr=paste0(wbChr[RefChr2],ifelse(RefStrand2=="True","+","-"),
                          ":",RefStart2,"-",RefEnd2))

# when seqContacts are interchromosomal, the stand, size and orientation of the hop are meaningless
seqContacts$hopSize[seqContacts$hopType=="trans"]<-NA
seqContacts$sameStrand[seqContacts$hopType=="trans"]<-NA
seqContacts$pairOri[seqContacts$hopType!="cis"]<-NA
seqContacts$pairOri[is.na(seqContacts$hopType)]<-NA

# correct hop orientation when fragments overlap
df<-seqContacts %>%
  dplyr::select(overlap,sameStrand,RefStrand1,RefStart1,RefStart2,pairOri,ReadID)
seqContacts$pairOri<-fixOverlapOri(df)

# remove any reads that do not have any seqContacts (not necessary?)
#readsWithseqContacts<-seqContacts$ReadID[!is.na(seqContacts$hopType)]
#idx<-seqContacts$ReadID %in% readsWithseqContactss
#hopReads<-seqContacts[idx,]

# keep only the hop information
#seqContactsOnly<-hopReads[!is.na(hopReads$hopType),]
# save table

importantStats["numContacts"]<-dim(seqContacts)[1]
importantStats["numCis"]<-sum(seqContacts$hopType=="cis")
importantStats["numTrans"]<-sum(seqContacts$hopType=="trans")
importantStats["numOverlap"]<-sum(seqContacts$overlap==TRUE)

write.csv(seqContacts,paste0("./csv/",sampleName,"_seqContacts.csv"),row.names=F)
saveRDS(seqContacts,paste0("./rds/",sampleName,"_seqContacts.rds"))


p1<-ggplot(seqContacts,aes(x=hopType,fill=hopType)) + geom_bar() + scale_fill_brewer(palette="Dark2") +
  ggtitle("Inter- vs intra-chromosomal (sequential contacts)")





p2<-seqContacts %>%
  filter(!is.na(pairOri)) %>%
  ggplot(aes(x=pairOri,fill=pairOri))  + geom_bar() + scale_fill_brewer(palette="Dark2") +
  ggtitle("Orientation of fragments (sequential contacts)")

p3<-seqContacts %>%
  filter(!is.na(pairOri)) %>%
  ggplot(aes(x=overlap,fill=overlap))  + geom_bar() + facet_wrap(~pairOri) +
  scale_fill_brewer(palette="Dark2") +
  ggtitle("Overlap of fragments (sequential contacts)")

p4<- seqContacts %>%
  filter(!is.na(pairOri)) %>%
  ggplot(aes(x=sameStrand,fill=sameStrand))  + geom_bar() + facet_wrap(~pairOri) +
  scale_fill_brewer(palette="Dark2") +
  ggtitle("Strandedness of fragments (sequential contacts)")

p5<- seqContacts %>%
  filter(!is.na(pairOri)) %>%
  ggplot(aes(x=sameStrand,y=hopSize,fill=sameStrand))  + geom_boxplot(outlier.shape = NA) +
  scale_y_continuous(limits = quantile(seqContacts$hopSize,na.rm=T, c(0.1, 0.9))) +
  scale_fill_brewer(palette="Dark2") +
  ggtitle("Distance between fragments (sequential contacts)")

p6<- seqContacts %>%
  filter(!is.na(pairOri)) %>%
  ggplot(aes(x=pairOri,y=hopSize,fill=pairOri))  + geom_boxplot(outlier.shape = NA) +
  scale_y_continuous(limits = quantile(seqContacts$hopSize,na.rm=T, c(0.1, 0.9))) +
  scale_fill_brewer(palette="Dark2") +
  ggtitle("Distance between fragments (sequential contacts)")

p7 <- seqContacts %>%
  filter(!is.na(pairOri)) %>%
  ggplot(aes(x=hopSize))  + geom_histogram(bins=60) +
  ggtitle("Distance between fragments (sequential contacts)")

minHop=500
p8 <- seqContacts %>%
  filter(!is.na(pairOri)) %>%
  filter(hopSize>minHop) %>%
  ggplot(aes(x=hopSize))  + geom_histogram(bins=60) +
  ggtitle(paste0("Distance between fragments (>",minHop,"bp), (sequential contacts)"))

p9 <- seqContacts %>%
  filter(!is.na(pairOri)) %>%
  ggplot(aes(x=hopSize,fill=sameStrand))  + geom_histogram(bins=60) + facet_wrap(~sameStrand) +
  scale_fill_brewer(palette="Dark2") +
  ggtitle(paste0("Distance between fragments by strandedness (sequential contacts)"))

p10 <- seqContacts %>%
  filter(!is.na(pairOri)) %>%
  ggplot(aes(x=hopSize,fill=pairOri))  + geom_histogram(bins=60) + facet_wrap(~pairOri) +
  scale_fill_brewer(palette="Dark2") +
  ggtitle(paste0("Distance between fragments by orientation (sequential contacts)"))


ml<-marrangeGrob(grobs=list(p1,p2,p3,p4,p5,p6,p7,p8,p9,p10),ncol=2,nrow=1)
ggsave(paste0("./plots_js/", sampleName,"_basicSeqContactsStats.pdf"),ml,device="pdf",
       width=29,height=20, units="cm")





######### look at data at the per read level

seqContacts<-seqContacts %>%
  group_by(ReadID) %>%
  mutate(numContacts=n())

# compare number of inter vs intra chromosomal contacts by number of contacts in a read
dt<-seqContacts %>%
  group_by(ReadID,numContacts,hopType) %>%
  summarise(count=n())

p1<-ggplot(dt,aes(x=as.factor(numContacts),y=count,fill=as.factor(hopType))) +
  geom_boxplot()  +
  ggtitle("Proportion of interchromosomal hops vs total number of hops") +
  scale_x_discrete(name="Contacts per read") + guides(fill=guide_legend(title="Hop type"))

# compare number of inter vs intra chromosomal contacts by read length
dt<-seqContacts %>%
  mutate(ReadLenRange=cut_width(ReadLength,width=1000,boundary=0)) %>%
  group_by(ReadLenRange,ReadID,hopType) %>%
  summarise(countHopType=n())

dt$ReadLenRange<-as.numeric(lapply(strsplit(as.character(dt$ReadLenRange),"[],]"),"[[",2))/1000

p2<-ggplot(dt,aes(x=as.factor(ReadLenRange),y=countHopType,fill=as.factor(hopType))) +
  geom_boxplot()  +
  ggtitle("Number of interchromosomal hops per read vs read length") +
  scale_x_discrete(name="Read length (kb)") + guides(fill=guide_legend(title="Hop type"))
p2
# compare number of inter vs intra chromosomal contacts per read
dt<-seqContacts %>%
  group_by(ReadID,ReadLength,hopType) %>%
  summarise(count=n()) %>%
  spread(key=hopType,value=count,fill=0)

p3<-ggplot(dt,aes(x=cis,y=trans)) +
  geom_jitter(col="darkblue",alpha=0.3) +
  ggtitle("Inter vs Intra sequential contacts per read") +
  scale_x_continuous(name="Number of Intra-chromosomal sequential contacts per read") +
  scale_y_continuous(name="Number of Inter-chromosomal sequential contacts per read")

ml<-marrangeGrob(grobs=list(p1,p2,p3),ncol=1,nrow=1)
ggsave(paste0(path, expName,"_sequentialContactsByRead.pdf"),ml,device="pdf",width=20,height=20,
       units="cm")

saveRDS(importantStats,paste0("./rds/importantStats_",sampleName,".rds"))

