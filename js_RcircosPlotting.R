#' Prepare reads for plotting with Rcircos
#'
#' Input data frame should have column names ReadID, Chr, RefStart and RefEnd
#' @param readData - data frame of fragments detected in MC-HiC
#' @param readsToDraw - vector of read IDs
#' @return data frame with 6 columns for pairs of interacting loci
#' @export
prepareLinkData<-function(readData,readsToDraw) {
  firstChr<-c()
  firstStart<-c()
  firstEnd<-c()
  secondChr<-c()
  secondStart<-c()
  secondEnd<-c()

  for (rd in readsToDraw) {
    currentRead<-readData[readData$ReadID==rd,]
    firstChr <- c(firstChr, currentRead$Chr[1:(dim(currentRead)[1]-1)])
    firstStart <- c(firstStart, currentRead$RefStart[1:(dim(currentRead)[1]-1)])
    firstEnd <- c(firstEnd, currentRead$RefEnd[1:(dim(currentRead)[1]-1)])
    secondChr <- c(secondChr, currentRead$Chr[2:(dim(currentRead)[1])])
    secondStart <- c(secondStart, currentRead$RefStart[2:(dim(currentRead)[1])])
    secondEnd <- c(secondEnd, currentRead$RefEnd[2:(dim(currentRead)[1])])
  }

  RCircosLink<- data.frame(firstChr=firstChr,firstStart=firstStart,firstEnd=firstEnd,
                           secondChr=secondChr,secondStart=secondStart,secondEnd=secondEnd,
                           stringsAsFactors=F)
  return(RCircosLink)
}

#' Prepare the core Rcircos plot for C. elegans data
#'
#' @param base.per.unit - integer for the size of the units that are plotted
#' @param chr.exclude - vector of names of chromosomes to exclude
#' @param track.inside - number of tracks to have inside the circle
#' @param track.outside - number of tracks to have outside the circle
#' @return plots ideogram
#' @export
baseRcircosCE<-function(base.per.unit=3000, chr.exclude=NULL, highlight.width=10, tracks.inside=1, tracks.outside=0){
  Chrnames<-c("chrI","chrII","chrIII","chrIV","chrV","chrX","MtDNA") # used to get rid of mtDNA
  ce11 <- list( "chrI" = 15072434,
                "chrII" = 15279421,
                "chrIII" = 13783801,
                "chrIV" = 17493829,
                "chrV" = 20924180,
                "chrX" = 17718942,
                "MtDNA" = 13794)

  ce11.ideo<-data.frame(Choromsome=Chrnames,ChromStart=0,ChromEnd=unlist(ce11),Band=1,Stain="gvar")
  cyto.info <- ce11.ideo
  RCircos.Set.Core.Components(cyto.info, chr.exclude,tracks.inside, tracks.outside)

  rcircos.params <- RCircos.Get.Plot.Parameters()
  rcircos.params$base.per.unit<-base.per.unit
  rcircos.params$chrom.width=0 #0.1
  rcircos.params$highlight.width=highlight.width #1
  RCircos.Reset.Plot.Parameters(rcircos.params)
  RCircos.Set.Plot.Area()
  par(mai=c(0.25, 0.25, 0.25, 0.25))
  plot.window(c(-1.5,1.5), c(-1.5, 1.5))

  RCircos.Chromosome.Ideogram.Plot()
}


#' Prepare a list of points of view for 4C
#'
#' Will use chromosome length to find positions at 20%, 50% and 80% of chromosome's
#' length to act as points of view for arms and center
#' @param chrLengthList - a named list with lengths of chromsomes
#' @param winSize - the size of the window around the POV for selecting interactions (must be an even number)
#' @return data.frame with points of view
#' @export
generatePOV<-function(chrLengthList=NULL,winSize=10000){
  if (is.null(chrLengthList)){
    chrLengthList <- list( "chrI" = 15072434,
                  "chrII" = 15279421,
                  "chrIII" = 13783801,
                  "chrIV" = 17493829,
                  "chrV" = 20924180,
                  "chrX" = 17718942)
    chrLengthList<-(unlist(chrLengthList))
  }
  left<-round(0.2*chrLengthList/1000,0)*1000
  center<-round(0.5*chrLengthList/1000,0)*1000
  right<-round(0.8*chrLengthList/1000,0)*1000
  names(left)<-paste(names(left),"left",sep="_")
  names(right)<-paste(names(right),"right",sep="_")
  names(center)<-paste(names(center),"center",sep="_")
  POV<-data.frame(POVname=c(names(left),names(center),names(right)),
                     POVpos=c(left,center,right),row.names=NULL)
  POV$chr<-gsub("_.*","",POV$POVname)
  POV$start<-POV$POVpos-winSize/2
  POV$end<-POV$POVpos+winSize/2
  POV<-POV[order(POV$chr,POV$start),]
  return(POV)
}
