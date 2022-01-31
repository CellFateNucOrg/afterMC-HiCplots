library(rhdf5)
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
chrOrder<-c("I","II","III","IV","V","X","MtDNA") # used to reorder sequence name levels in objects
importantStats<-list()
dir.create("./plots_js")
dir.create("./rds")
dir.create("./csv")


sampleName="20190501_HIC6_7_barcode08_pass"
fileName=paste0("./frag_files/frg_",sampleName,"_WS235.hdf5")
subSample="raw"
####################
#### read in data
####################
h5ls(fileName)

myData<- as_tibble(t(h5read(fileName,"frg_np")),.name_repair="minimal")
names(myData)<-h5read(fileName,"frg_np_header_lst")
dim(myData)

####################
#### general stats
####################
importantStats["TotalReadsSequenced"]<-max(myData$ReadID)
importantStats["TotalFragsSequenced"]<-max(myData$FrgID)
importantStats["TotalMappedReads"]<-length(unique(myData$ReadID))
myData$uniqFragID<-paste(myData$ReadID,myData$FrgID,sep="_")
importantStats["TotalMappedFrags"]<-dim(myData)[1]

####################
#### filtering reads
####################
multiMapFragID<-unique(myData$uniqFragID[duplicated(myData$uniqFragID)])
singleMapped<-myData[-which(myData$uniqFragID %in% multiMapFragID),]
singleMapped$uniqueMapping<-1
importantStats["uniquelyMappedFrags"]<-dim(singleMapped)[1]

# filter multimapped fragments by choosing best MQ, or taking first in list
# uniqueMapping score is 1 if unique, 0 if multimapped but has single frag
# with highest MQ,
multiMappedFrags<-myData[which(myData$uniqFragID %in% multiMapFragID),]
multiMappedFrags$uniqueMapping<-0
importantStats["multiMappedFrags"]<-dim(multiMappedFrags)[1]


# take alignment with best MQ
mm<-multiMappedFrags %>% group_by(uniqFragID) %>% filter(MQ==max(MQ))
# remove ambiguous alignments where more than one has top MQ
ambig<-unique(mm[duplicated(mm$uniqFragID),]$uniqFragID)
importantStats["ambigMappedFrags"]<-length(ambig)

mm<-mm[! mm$uniqFragID %in% ambig,]

singleMapped<-bind_rows(singleMapped,mm)
# re-sort tibble
singleMapped<-singleMapped[order(singleMapped$ReadID,singleMapped$FrgID),]
singleMapped$FragLength<-abs(singleMapped$SeqStart-singleMapped$SeqEnd)

# separate mtDNA frags
withMtDNA<-singleMapped
### save same data with mtDNA (true hop and MQ>=20 only)
saveRDS(withMtDNA,paste0("./rds/withMtDNA_",sampleName,".rds"))

# analyse data with out mtDNA
singleMapped<-singleMapped[!singleMapped$Chr==7,]



importantStats["bestUniqMappedFrags"]<-dim(singleMapped)[1]
importantStats["trueHops"]<-dim(singleMapped[singleMapped$TrueHop==1,])[1]
importantStats["trueHops_MQgt20"]<-dim(singleMapped[singleMapped$TrueHop==1 & singleMapped$MQ>=20,])[1]
importantStats["medianReadLength"]<-median(singleMapped$ReadLength,na.rm=T)
importantStats["medianFragLength"]<-median(singleMapped$FragLength,na.rm=T)
importantStats["minReadLength"]<-min(singleMapped$ReadLength,na.rm=T)
importantStats["minFragLength"]<-min(singleMapped$FragLength,na.rm=T)
importantStats["maxReadLength"]<-max(singleMapped$ReadLength,na.rm=T)
importantStats["maxFragLength"]<-max(singleMapped$FragLength,na.rm=T)



####################
#### MQ plots
####################
p0<-ggplot(singleMapped, aes(x=MQ)) + geom_histogram(bins=125) +
  geom_vline(xintercept=20, col="red") +
  ggtitle("Histogram of mapping quality for uniquely & best mapping reads")

p1<-ggplot(singleMapped,aes(x=MQ,y=ReadLength)) + geom_hex(bins=100) +
  geom_vline(xintercept=20, col="red") +
  ggtitle("Mapping quality vs read length") +
  scale_fill_viridis()

p2<-ggplot(singleMapped,aes(x=MQ,y=FragLength)) + geom_hex(bins=100) +
  geom_vline(xintercept=20, col="red") +
  ggtitle("Mapping quality vs fragment length") +
  scale_fill_viridis()


p<-ggarrange(plotlist=list(p0,p1,p2),nrow=3)
ggsave(paste0("./plots_js/mappingQuality_",sampleName,"_",subSample,".pdf"),plot=p, device="pdf", width=17, height=26,units="cm")




####################
#### fragments per read
####################





### all reads

subSample="bestUniqMapped"
readData <- singleMapped %>%
  group_by(ReadID) %>%
  summarise(
    numFrags = length(unique(FrgID)),
    numChr = length(unique(Chr)),
    readLength = unique(ReadLength))

importantStats[paste0("maxFragPerRead_",subSample)]<-max(readData$numFrags)
importantStats[paste0("FractionMultiFragReads_",subSample)]<-round(sum(readData$numFrags>1)/dim(readData)[1],2)

p<-marrangeGrob(plotPerReadData(readData,subSample),nrow=2,ncol=1)

ggsave(paste0("./plots_js/plotsPerRead_",sampleName,"_",subSample,".pdf"),
       plot=p, device="pdf",
       width=17,height=26,units="cm")



### reads with MQ>=20
subSample="MQ>=20"
readData <- singleMapped[singleMapped$MQ>=20,] %>%
  group_by(ReadID) %>%
  summarise(
    numFrags = length(unique(FrgID)),
    numChr = length(unique(Chr)),
    readLength = unique(ReadLength))

importantStats[paste0("maxFragPerRead_",subSample)]<-max(readData$numFrags)
importantStats[paste0("FractionMultiFragReads_",subSample)]<-round(sum(readData$numFrags>1)/dim(readData)[1],2)

p<-marrangeGrob(plotPerReadData(readData,subSample),nrow=2,ncol=1)

ggsave(paste0("./plots_js/plotsPerRead_",sampleName,"_",subSample,".pdf"),
       plot=p, device="pdf",width=17,height=26,units="cm")


### reads with TrueHop=TRUE
subSample="TrueHop"
readData <- singleMapped[singleMapped$TrueHop==1 & singleMapped$MQ>=20,] %>%
  group_by(ReadID) %>%
  summarise(
    numFrags = length(unique(FrgID)),
    numChr = length(unique(Chr)),
    readLength = unique(ReadLength))

importantStats[paste0("maxFragPerRead_",subSample)]<-max(readData$numFrags)
importantStats[paste0("FractionMultiFragReads_",subSample)]<-round(sum(readData$numFrags>1)/dim(readData)[1],2)

p<-marrangeGrob(plotPerReadData(readData,subSample),nrow=2,ncol=1)


tbl<-table(readData$numFrags[readData$numFrags>1])
sum(tbl)/sum(readData$numFrags)
tbl<-table(readData$numFrags[readData$numFrags>2])
sum(tbl)/sum(readData$numFrags)

ggsave(paste0("./plots_js/plotsPerRead_",sampleName,"_",subSample,".pdf"),
       plot=p, device="pdf",width=17,height=26,units="cm")

#### save stats
saveRDS(importantStats,paste0("./rds/importantStats_",sampleName,".rds"))

singleMapped <- singleMapped[singleMapped$TrueHop==1 & singleMapped$MQ>=20,] %>%
  group_by(ReadID) %>%
  mutate(numFrags = length(unique(FrgID)),
         numChrs = length(unique(Chr)))

# withMtDNA <- withMtDNA[withMtDNA$TrueHop==1 & withMtDNA$MQ>=20,] %>%
#   group_by(ReadID) %>%
#   mutate(numFrags = length(unique(FrgID)),
#          numChrs = length(unique(Chr)))

### save bestUniqMapped (true hop and MQ>=20 only)
saveRDS(singleMapped,paste0("./rds/bestUniqMapped_",sampleName,".rds"))

### save same data with mtDNA (true hop and MQ>=20 only)
#saveRDS(withMtDNA,paste0("./rds/withMtDNAbest_",sampleName,".rds"))


### save fragments from multiFrag reads
saveRDS(singleMapped[singleMapped$numFrags>1,],
        paste0("./rds/multiFrag_",sampleName,".rds"))
write.table(singleMapped[singleMapped$numFrags>1,],
        paste0("./csv/multiFrag_",sampleName,".tsv"),
        sep="\t",quote=F, row.names=F)




saveRDS(importantStats,paste0("./rds/importantStats_",sampleName,".rds"))

########
# histogram of read length by Fragment number
########
subSample="best"
bestUniqMapped<-readRDS(paste0("./rds/bestUniqMapped_",sampleName,".rds"))

readData<- bestUniqMapped %>%
  distinct(ReadID,.keep_all=T) %>%
  filter(numFrags<5) %>%
  group_by(numFrags)

p0<-ggplot(readData, aes(x=ReadLength)) + geom_histogram(bins=125) +
  geom_vline(xintercept=c(500,1000,2000), col="red") +
  ggtitle("Histogram of read length") + facet_wrap(~numFrags) +xlim(0,5000)

ggsave(paste0("./plots_js/readLengthByFragNum_",sampleName,"_",subSample,".pdf"),plot=p0, device="pdf", width=17, height=17,units="cm")

readData<- bestUniqMapped %>%
  distinct(ReadID,.keep_all=T) %>%
  mutate(bin=ifelse(floor(ReadLength/1000)>5,5,floor(ReadLength/1000))) %>%
  mutate(numFragsBin=ifelse(numFrags>8,8,numFrags))

p1<-ggplot(readData) + geom_bar(aes(x=as.factor(bin),fill=as.factor(numFragsBin)),position="stack") +
  scale_fill_brewer(palette="Dark2") + labs(fill="Fragments\nper read",y="Number of reads") +
  ggtitle("Number of fragments in reads of different lengths")  +
  scale_x_discrete(name="Read length (kb)", breaks=c(0,1,2,3,4,5),labels=c("0-1","1-2","2-3","3-4","4-5",">6"))

p2<-ggplot(readData) + geom_bar(aes(x=as.factor(bin),fill=as.factor(numFragsBin)),position="fill") +
  scale_fill_brewer(palette="Dark2") + labs(fill="Fragments\nper read",y="Fraction of reads") +
  ggtitle("Number of fragments in reads of different lengths")  +
  scale_x_discrete(name="Read length (kb)", breaks=c(0,1,2,3,4,5),labels=c("0-1","1-2","2-3","3-4","4-5",">6"))

p<-ggarrange(plotlist=list(p1,p2),nrow=2)
ggsave(paste0("./plots_js/fragmentsPerRead_",sampleName,"_",subSample,".pdf"),plot=p, device="pdf", width=17, height=26,units="cm")




########
# include mitochondrial data as control
########
subSample="MtDNA"
withMtDNA<-readRDS(paste0("./rds/withMtDNA_",sampleName,".rds"))

importantStats["numFragsMtDNA"]<-sum(withMtDNA$Chr==7)
importantStats["percentFragsMtDNA"]<-100*sum(withMtDNA$Chr==7)/dim(withMtDNA)[1]
# Note: MtDNA size is 13794 bp
importantStats["percentGenomeMtDNA"]<-100*13794/1e8

readData<-withMtDNA %>%
  group_by(ReadID) %>%
  summarise(numFrags = length(unique(FrgID)),
            mtDNAcount = sum(Chr==7),
            nonMtDNAcount = sum(Chr!=7)) %>%
  filter(mtDNAcount>0)

# number of reads containing mtDNA fragment
importantStats["numAllReads"]<-length(unique(withMtDNA$ReadID))
importantStats["numReadsMtDNA"]<-dim(readData)[1]
importantStats["percentReadsWithMtDNA"]<-round(100*importantStats$numReadsMtDNA/importantStats$numAllReads,2)

readData<-withMtDNA %>%
  group_by(ReadID) %>%
  summarise(numFrags = length(unique(FrgID)),
            mtDNAcount = sum(Chr==7),
            nonMtDNAcount = sum(Chr!=7)) %>%
  filter(mtDNAcount>0) %>%
  filter(nonMtDNAcount>0)

importantStats["numMixedReadsMtDNA"]<-dim(readData)[1]
importantStats["percentMixedReadsWithMtDNA"]<-round(100*importantStats$numMixedReadsMtDNA/importantStats$numAllReads,2)

saveRDS(importantStats,paste0("./rds/importantStats_",sampleName,".rds"))



#############################################################
#### RCircos
############################################################

#library(circlize)
library("RCircos")
library(BSgenome.Celegans.UCSC.ce11)
library(dplyr)
library(tidyr)

source("./js_RcircosPlotting.R")


subSample="MtDNA"
withMtDNA<-readRDS(paste0("./rds/withMtDNA_",sampleName,".rds"))

readData<-withMtDNA %>%
  group_by(ReadID) %>%
  mutate(mtDNAcount = sum(Chr==7)) %>%
  mutate(numFrags = n()) %>%
  filter(mtDNAcount>0) %>%
  filter(numFrags>1)

Chrnames<-c("chrI","chrII","chrIII","chrIV","chrV","chrX","MtDNA") # used to get rid of mtDNA

############
## Prepare reads to plot
############

readData$Chr<-Chrnames[readData$Chr]

readsToDraw<-unique(readData$ReadID)
#readsToDraw<-sample(unique(readData$ReadID),1)
#readsToDraw<-c(102)



RCircosLink<-prepareLinkData(readData,readsToDraw)

pdf(file=paste0("./plots_js/circos_",sampleName,"_",subSample,".pdf"),width=8,height=8,
    paper="a4",title="Reads with mtDNA fragments")

baseRcircosCE()
track.num <- 1
RCircos.Link.Plot(RCircosLink, track.num=track.num, by.chromosome=TRUE,is.sorted=FALSE)
dev.off()





#RCircos.Ribbon.Plot(ribbon.data=RCircos.Ribbon.Data,
#                    + track.num=11, by.chromosome=FALSE, twist=FALSE)

