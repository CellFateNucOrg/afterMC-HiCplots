# Using domain definitions from:
#"Caenorhabditis elegans chromosome arms are anchored to the nuclear membrane via discontinuous association with LEM-2"
#Kohta Ikegami, Thea A Egelhofer, Susan Strome & Jason D Lieb
#Genome Biology volume 11, Article number: R120 (2010)
# https://genomebiology.biomedcentral.com/articles/10.1186/gb-2010-11-12-r120
# Table S7 in additional file 2

options(scipen=9)
domains<-read.delim("chrDomains.tsv",header=T,stringsAsFactors=F)
binSize=50000
#round to bin size

#ce11<-read.delim("ce11.chrom.sizes",header=F,stringsAsFactors=F)
#names(ce11)<-c("chr","size")

#idx<-match(domains$chromosome,ce11$chr)
#domains$ce11size<-ce11$size[idx]
#domains$biggestSize<-domains$end
#domains$biggestSize[domains$domain=="Right"]<-apply(cbind(
#  domains$end[domains$domain=="Right"],
#  domains$ce11[domains$domain=="Right"]),1,max)

domains$domainName<-paste(domains$chromosome, domains$domain,sep="_")
#domains$pseudoStart<-domains$start-domains$start
#domains$pseudoEnd<-domains$end-domains$start+1
domains$roundedStart<-ceiling((domains$start-1)/binSize)*binSize
domains$roundedEnd<-ceiling(domains$end/binSize)*binSize

domains$pseudoStart<-0
domains$pseudoEnd<-(ceiling((domains$roundedEnd-domains$roundedStart)/binSize))*binSize



# convert to ordinal numbers for todor numbering
domains$chr<-as.factor(domains$chromosome)
levels(domains$chr)<-1:6
domains$ranges<-paste0(domains$chr,":",
                       apply(cbind(domains$roundedStart,domains$start),1,max),"-",
                       apply(cbind(domains$roundedEnd, domains$end),1,min))


chrom.sizes<-domains[domains$domain=="Left",]
write.table(cbind(chrom.sizes$domainName, chrom.sizes$pseudoEnd),
            file="left.chrom.sizes",quote=F,sep="\t",row.names=F,col.names=F)
write.table(cbind(chrom.sizes$ranges, chrom.sizes$domainName, chrom.sizes$roundedStart),
            file="left.chrom.ranges", quote=F, sep="\t", row.names=F, col.names=F)


chrom.sizes<-domains[domains$domain=="Right",]
write.table(cbind(chrom.sizes$domainName, chrom.sizes$pseudoEnd),
            file="right.chrom.sizes",quote=F,sep="\t",row.names=F,col.names=F)
write.table(cbind(chrom.sizes$ranges, chrom.sizes$domainName, chrom.sizes$roundedStart),
            file="right.chrom.ranges", quote=F, sep="\t", row.names=F, col.names=F)

chrom.sizes<-domains[domains$domain=="Center",]
write.table(cbind(chrom.sizes$domainName, chrom.sizes$pseudoEnd),
            file="center.chrom.sizes",quote=F,sep="\t", row.names=F, col.names=F)
write.table(cbind(chrom.sizes$ranges, chrom.sizes$domainName, chrom.sizes$roundedStart),
            file="center.chrom.ranges", quote=F, sep="\t", row.names=F, col.names=F)

