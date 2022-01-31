#https://cran.r-project.org/web/packages/BioCircos/vignettes/BioCircos.html
#devtools::install_github('lvulliard/BioCircos.R', build_vignettes = TRUE)
library(BioCircos)
library(BSgenome.Celegans.UCSC.ce11)


sampleName="20190501_HIC6_7_barcode08_pass"
multiFrag<-readRDS(paste0("./rds/multiFrag_",sampleName,".rds"))
Chrnames<-c("chrI","chrII","chrIII","chrIV","chrV","chrX","MtDNA") # used to get rid of mtDNA
ce11 <- list( "chrI" = 15072434,
              "chrII" = 15279421,
              "chrIII" = 13783801,
              "chrIV" = 17493829,
              "chrV" = 20924180,
              "chrX" = 17718942,
              "MtDNA" = 13794)

# BioCircos(genome = ce11, genomeFillColor = c("tomato2", "darkblue"),
#           genomeTicksScale = 1e+6, genomeLabelTextSize=28,  genomeTicksLen = 2)
#
#
#
# gt4<-multiFrag[multiFrag$numFrags==4,]
#
# gt4$Chr<-Chrnames[gt4$Chr]
#
# seed(2020)
# readsToDraw<-sample(unique(gt4$ReadID),8)
#
# plotList<-list()
#
# for (rd in readsToDraw) {
#   currentRead<-gt4[gt4$ReadID==rd,]
#   links_firstChr<-currentRead$Chr[1:3]
#   links_secondChr<-currentRead$Chr[2:4]
#   links_pos_firstChr<-currentRead$RefStart[1:3]
#   links_pos_secondChr<-currentRead$RefEnd[2:4]
#   tracklist = BioCircosBackgroundTrack("myBackgroundTrack", minRadius = 0, maxRadius = 0.9,
#                                        borderSize = 0, fillColors = "#EEFFEE")
#
#   tracklist = tracklist + BioCircosLinkTrack('myLinkTrack',
#                                              links_firstChr,links_pos_firstChr,
#                                              links_pos_firstChr,
#                                              links_secondChr, links_pos_secondChr,
#                                               links_pos_secondChr,
#                                              maxRadius = 0.9)
#
#   BioCircos(genome=ce11, tracklist, genomeFillColor = "PuOr",
#             chrPad = 0.02, displayGenomeBorder = FALSE,
#             genomeTicksDisplay = FALSE,  genomeLabelTextSize = "22pt", genomeLabelDy = 0)
# }
#
# rd=117




#############################################################
#### RCircos
############################################################

#####################
## plotting a point of view
#####################


#library(circlize)
library("RCircos")
library(BSgenome.Celegans.UCSC.ce11)

source("./js_RcircosPlotting.R")

sampleName="20190501_HIC6_7_barcode08_pass"
multiFrag<-readRDS(paste0("./rds/multiFrag_",sampleName,".rds"))

Chrnames<-c("chrI","chrII","chrIII","chrIV","chrV","chrX","MtDNA") # used to get rid of mtDNA


winSizes=c(5000,10000,20000,50000)

for (winSize in winSizes){
  pdf(file=paste0("./plots_js/circos_",sampleName,"_",subSample,"_win",winSize/1000,"kb.pdf"),
      width=8, height=11, paper="a4")
  par(mfrow=c(4,3))

  POV<-generatePOV(winSize=winSize)
  POV$chr<-as.integer(factor(POV$chr))

  for (pov in 1:dim(POV)[1]) {
    subSample<-"POV"
    readData<-multiFrag %>%
      filter(Chr==POV$chr[pov]) %>%
      filter(RefStart>POV$start[pov]) %>%
      filter(RefEnd<POV$end[pov])

    readsToDraw<-unique(readData$ReadID)

    readData<-multiFrag %>%
      filter(ReadID %in% readsToDraw)

    readData$Chr<-Chrnames[readData$Chr]

    #RCircosLink<-prepareLinkData(readData,readsToDraw)
    baseRcircosCE(chr.exclude = "MtDNA",highlight.width=3)
    for (rd in readsToDraw){
      RCircosLink<-prepareLinkData(readData,rd)
      track.num <- 1
      RCircos.Link.Plot(RCircosLink, track.num, by.chromosome=TRUE)
    }
  }
  dev.off()
}



#####################
## plotting a sample of single reads
#####################


#library(circlize)
library("RCircos")
library(BSgenome.Celegans.UCSC.ce11)

source("./js_RcircosPlotting.R")

sampleName="20190501_HIC6_7_barcode08_pass"
multiFrag<-readRDS(paste0("./rds/multiFrag_",sampleName,".rds"))

Chrnames<-c("chrI","chrII","chrIII","chrIV","chrV","chrX","MtDNA") # used to get rid of mtDNA

for (fragNum in 3:8) {
  subSample<-paste0(fragNum,"fragReads")
  readData<-multiFrag[multiFrag$numFrags==fragNum,]
  readData$Chr<-Chrnames[readData$Chr]

  #readsToDraw<-unique(readData$ReadID)
  set.seed(2020)
  readsToDraw<-sample(unique(readData$ReadID),36)
  plottedReads[[paste0(fragNum,"fragReads")]]<-readsToDraw
  #readsToDraw<-c(102)

  #RCircosLink<-prepareLinkData(readData,readsToDraw)

  pdf(file=paste0("./plots_js/circos_",sampleName,"_",subSample,".pdf"),width=8,height=11,
      paper="a4")
  par(mfrow=c(4,3))

  for (rd in readsToDraw){
    baseRcircosCE(chr.exclude = "MtDNA",highlight.width=3)
    RCircosLink<-prepareLinkData(readData,rd)
    track.num <- 1
    RCircos.Link.Plot(RCircosLink, track.num, by.chromosome=TRUE)
  }
  dev.off()
}

write.csv(as.data.frame(plottedReads),file=paste0("./csv/circos_",sampleName,"_readsPlotted.csv"))
