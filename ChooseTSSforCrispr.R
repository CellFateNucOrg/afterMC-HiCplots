library(BSgenome.Celegans.UCSC.ce11)
library(Biostrings)
library(GenomicRanges)

########################################################################################

##################
### functions
##################

miObj2gr<-function(miObj) {
  #function to convert Mindex object (obtained from applying GR on DNAstringset) to genomic ranges
  # this function is required by call_context_methylation
  allGR<-GRanges()
  seqlevels(allGR)<-names(miObj)
  for (n in names(miObj)) {
    if (length(miObj[[n]])>0) {
      grObj<-GRanges(seqnames=Rle(c(n),length(miObj[[n]])),
                     ranges=miObj[[n]], strand="*")
      allGR<-append(allGR,grObj)
    }
  }
  return(allGR)
}

########################################################################################







# read in the TSSs for amplicons
amptss<-readRDS("/Users/semple/Documents/MeisterLab/dSMF/PromoterPrimerDesign/usefulObjects/ampliconMaxTSSgr.RDS")

fragSize=600

# expand to region of +-fragSize/2 around TSS (must NOT contain DpnII site)
tss<-resize(amptss,width=fragSize,fix="center")

#restriction enzyme motifs
DpnII<-"GATC"
NlaIII<-"CATG"

# get sequence of Clegans genome
genome<-getSeq(Celegans)

# find dpnII sites
dpnSites<-miObj2gr(vmatchPattern(DpnII,genome))


# find NlaIII sites
nlaSites<-miObj2gr(vmatchPattern(NlaIII,genome))

# remove TSSs that have dpnII sites within +-fragSize/2  of them
ol<-findOverlaps(dpnSites,tss,ignore.strand=TRUE)
nonDpnSites<-tss[-unique(subjectHits(ol))]
length(nonDpnSites)


# remove TSSs that do not contain nlaIII sites within +-fragSize/2  of them
ol<-findOverlaps(nlaSites,nonDpnSites)
withNlaSites<-nonDpnSites[unique(subjectHits(ol))]
length(withNlaSites)


# export the sequences so you can feed them into benchling
seqs<-getSeq(Celegans, withNlaSites)
names(seqs)<-withNlaSites$WBgeneID
writeXStringSet(seqs, filepath=paste0("~/Downloads/fasta_tss",fragSize,".fa"), format="fasta")

