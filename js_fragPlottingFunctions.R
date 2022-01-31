#####################
##### plot per read functions
####################


plotPerReadData<-function(readData,subSample) {
  # plot number of fragments per read
  dt <- readData %>%
    group_by(numFrags) %>%
    summarise(Count = n()) %>%
    mutate(Frequency=Count/sum(Count))

  p1 <- ggplot(readData) +
    geom_bar(aes(x = as.factor(numFrags)), fill = "darkblue") +
    ggtitle(paste0("Number of fragments per read (", subSample, ")")) +
    scale_x_discrete(name = "Fragments per read",limits=1:15) +
    theme_classic()

  p2 <- ggplot(dt) +
    geom_col(aes(x = as.factor(numFrags),y=Frequency), fill = "darkgreen") +
    ggtitle(paste0("Number of fragments per read (", subSample, ")")) +
    scale_x_discrete(name = "Fragments per read",limits=1:15) +
    theme_classic()

  # plot number of fragments per read vs read length
  dt <- readData %>%
    mutate(readLength=cut_width(readLength/1000,width=2,boundary=0))

  p3<-ggplot(dt, aes(x = readLength, y = numFrags)) +
    geom_boxplot(fill = "lightblue",width=0.5) +
    theme(axis.text.x=element_blank()) +
    ggtitle(paste0("Number of fragments per read length (", subSample, ")")) +
    scale_x_discrete(name = "Read length (kb)") + theme_classic() +
    ylab("Number of fragments per read")

  dt <- readData %>%
    mutate(readLength=cut_width(readLength/1000,width=2,boundary=0))

  p4<-ggplot(dt, aes(x = readLength, y = numFrags)) +
    geom_boxplot(fill = "lightblue",varwidth=T) +
    theme(axis.text.x=element_blank()) +
    ggtitle(paste0("Number of fragments per read length (", subSample, ")")) +
    scale_x_discrete(name = "Read length (kb)") + theme_classic() +
    ylab("Number of fragments per read")

  # plot number of chromosomes per read
  dt <- readData %>%
    filter(numFrags>1) %>%
    group_by(numChr) %>%
    summarise(Count = n()) %>%
    mutate(Frequency=Count/sum(Count))

  p5<-ggplot(dt,aes(x=as.factor(numChr),y=Count)) +
    geom_bar(stat="identity",fill = "darkblue") +
    ggtitle(paste0("Chromosomes per read of multi-fragment reads (", subSample, ")")) +
    scale_x_discrete(name="Chromosomes per read") + theme_classic()

  p6<-ggplot(dt,aes(x=as.factor(numChr),y=Frequency)) +
    geom_bar(stat="identity",fill = "darkgreen") +
    ggtitle(paste0("Chromosomes per read of multi-fragment reads (", subSample, ")")) +
    scale_x_discrete(name="Chromosomes per read") + theme_classic()

  return(list(p1,p2,p3,p4,p5,p6))
}







#####################
##### consecutive hop functions
####################



minDistance<-function(v1,v2,v3,v4) {
  # function finds the smallest distance between all fragment ends
  #it takes four vectors representing AlnStart, AlnEnd, lead(AlnStart), lead(AlnEnd)
  # output is a vector with minimum absolute distance
  allDist<-cbind(v1-v3,v2-v4,v2-v3,v1-v4)
  minDist<-apply(abs(allDist),1,min)
  return(minDist)
}

pairOrientation<-function(df) {
  # function finds the relative orientation of two fragments
  # input is a dataframe with six columns: RefStart,RefEnd, lead(RefStart),lead(RefEnd) Strand, lead(Strand)
  # output is one of "divergent","convergent","tandem","inverted".
  df1<-df
  # reverse orientation of start end according to strand
  idx<-which(df[,5]=="-1")
  df1[idx,1]<-df[idx,2]
  df1[idx,2]<-df[idx,1]
  idx<-which(df[,6]=="-1")
  df1[idx,3]<-df[idx,4]
  df1[idx,4]<-df[idx,3]
  df<-df1
  allDist<-cbind(as.numeric(df[,1])-as.numeric(df[,3]),
                 as.numeric(df[,2])-as.numeric(df[,4]),
                 as.numeric(df[,2])-as.numeric(df[,3]),
                 as.numeric(df[,1])-as.numeric(df[,4]))
  # find rows with NA
  naRows<-rowSums(is.na(allDist))>0
  minType<-apply(abs(allDist[!naRows,]),1,which.min)
  ori<-c("divergent","convergent","tandem","inverted")[minType]
  allOri<-rep(NA,length(naRows))
  allOri[!naRows]<-ori
  return(allOri)
}


isOverlapping<-function(df) {
  # function finds fragment pairs from contacts that overlap
  # output is TRUE (overlapping) or FALSE (not overlapping)
  GR1<-GRanges(seqnames=df[,7],
               ranges=IRanges(start=as.numeric(df[,1]),
                              end=as.numeric(df[,2])), strand="*")
  GR2<-GRanges(seqnames=df[,8],
               ranges=IRanges(start=as.numeric(df[,3]),
                              end=as.numeric(df[,4])), strand="*")
  ol<-findOverlaps(GR1,GR2)
  idx<-which(queryHits(ol)==subjectHits(ol))
  overlapping<-rep(FALSE, length(GR1))
  overlapping[queryHits(ol)[idx]]<-TRUE
  return(overlapping)
}


getOverlapSize<-function(df) {
  # function finds fragment pairs from contacts that overlap
  # output is TRUE (overlapping) or FALSE (not overlapping)
  GR1<-GRanges(seqnames=df[,5],
               ranges=IRanges(start=as.numeric(df[,1]),
                              end=as.numeric(df[,2])), strand="*")
  GR2<-GRanges(seqnames=df[-dim(df)[1],6],
               ranges=IRanges(start=as.numeric(df[-dim(df)[1],3]),
                              end=as.numeric(df[-dim(df)[1],4])), strand="*")
  ol<-findOverlaps(GR1,GR2)
  idx<-which(queryHits(ol)==subjectHits(ol))
  sizeOverlap<-rep(0, length(GR1))
  sizeOverlap[queryHits(ol)[idx]]<-width(pintersect(GR1[queryHits(ol)[idx]],GR2[subjectHits(ol)[idx]]))
  return(sizeOverlap)
}



fixOverlapOri<-function(df) {
  # this function fixes the orientation calculated by pairOrientation which can fail if fragments
  # overlap.
  # output is vector with orientations
  names(df)<-c("overlap","sameStrand","RefStrand1","RefStart1","RefStart2","pairOri","ReadID")
  pairOri<-df$pairOri
  df<-df[,-6]
  i<-which(df$overlap==T & df$sameStrand==T & df$RefStrand1=="1")
  df1<-df[i,]
  df1$pairOri<-with(df1,
                    ifelse(RefStart1<RefStart2,"tandem","inverted"))
  pairOri[i]<-df1$pairOri
  ##
  i<-which(df$overlap==T & df$sameStrand==T & df$RefStrand1=="-1")
  df1<-df[i,]
  df1$pairOri<-with(df1,
                    ifelse(RefStart1>RefStart2,"tandem","inverted"))
  pairOri[i]<-df1$pairOri
  ##
  i<-which(df$overlap==T & df$sameStrand==F & df$RefStrand1=="1")
  df1<-df[i,]
  df1$pairOri<-with(df1,
                    ifelse(RefStart1<RefStart2,"convergent","divergent"))
  pairOri[i]<-df1$pairOri
  ##
  i<-which(df$overlap==T & df$sameStrand==F & df$RefStrand1=="-1")
  df1<-df[i,]
  df1$pairOri<-with(df1,
                    ifelse(RefStart1>RefStart2,"convergent","divergent"))
  pairOri[i]<-df1$pairOri
  ##
  return(pairOri)
}



areConsecutive<-function(FragmentId,FragmentId1){
  # checks if two fragments ids are consecutive, accepting complex ids formed
  # by merging several fragments (e.g 2.4 means merging of fragments 2-3-4)
  first=FragmentId%%1
  compound<-first!=0
  first[!compound]<-(FragmentId%/%1)[!compound]
  first[compound]=first[compound]*10
  second=FragmentId1%/%1
  return(second-first==1)
}
