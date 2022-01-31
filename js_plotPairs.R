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
#library(InteractionSet)
library(readr)
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

sampleName="SRR1665089_wt_N2_Crane2015"
fileName=paste0("./pairs/",sampleName,".bsorted.pairs")
subSample="raw"
####################
#### read in data
####################
#system(paste0("gunzip ",fileName,".gz"))

myData<-read_tsv(fileName, skip=12, col_names=FALSE,
                 col_types=list(col_character(),
                                col_character(),
                                col_integer(),
                                col_character(),
                                col_integer(),
                                col_character(),
                                col_character()))
names(myData)<-c("readID","chrom1","start1","chrom2","start2","strand1","strand2")
dim(myData)

####################
#### general stats
####################

importantStats["TotalMappedReads"]<-dim(myData)[1]

####################
#### filtering reads
####################
multiMapReadID<-unique(myData$readID[duplicated(myData$readID)])
singleMapped<-myData[-which(myData$readID %in% multiMapReadID),]
importantStats["uniquelyMappedReads"]<-dim(singleMapped)[1]


singleMapped$contactLength<-abs(singleMapped$start1-singleMapped$start2)
singleMapped$sameChr<-ifelse(singleMapped$chrom1==singleMapped$chrom2,"cis","trans")

# separate mtDNA frags
withMtDNA<-singleMapped
### save same data with mtDNA
saveRDS(withMtDNA,paste0("./rds/withMtDNA_",sampleName,".rds"))
singleMapped<-readRDS(paste0("./rds/withMtDNA_",sampleName,".rds"))

# analyse data with out mtDNA
singleMapped<-singleMapped[!(singleMapped$chrom1=="chrM" | singleMapped$chrom2=="chrM"),]

importantStats["numNuclearReads"]<-dim(singleMapped)[1]
importantStats["medianContactLength"]<-median(singleMapped$contactLength,na.rm=T)
importantStats["minContactLength"]<-min(singleMapped$contactLength,na.rm=T)
importantStats["maxContactLength"]<-max(singleMapped$contactLength,na.rm=T)
importantStats["numCisContacts"]<-sum(singleMapped$sameChr=="cis")
importantStats["numTransContacts"]<-sum(singleMapped$sameChr=="trans")
importantStats["fractionCisContacts"]<-importantStats[["numCisContacts"]]/dim(singleMapped)[1]

importantStats["numHopGT500"]<-length(singleMapped$contactLength>500




#### save stats
saveRDS(importantStats,paste0("./rds/importantStats_",sampleName,".rds"))


########
# histogram of contact length
########


p0<-ggplot(singleMapped, aes(x=log(contactLength))) + geom_histogram(bins=125) +
  geom_vline(xintercept=c(500,1000,2000), col="red") +
  ggtitle("Histogram of log contact length") #+xlim(0,5000)

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

