### 27 Sep 2016
### Katarzyna Wreczycka
### Some examples of how to use genomation


# library(devtools)
# install_github("BIMSBbioinfo/genomation",build_vignettes=FALSE)
library(genomation)
# source("http://bioconductor.org/biocLite.R")
# biocLite("genomationData")
library(genomationData)

# Read targets
genomationDataPath = system.file('extdata',package='genomationData')
bam.files = list.files(genomationDataPath, full.names=TRUE, pattern='bam$')
bam.files = bam.files[!grepl('Cage', bam.files)]
bam.files = bam.files[1:4]
bam.files
# Read windows
ctcf.peaks = readGeneric(file.path(genomationDataPath, 
                          'wgEncodeBroadHistoneH1hescCtcfStdPk.broadPeak.gz'),
                         meta.cols=list(name=4,
                                        score=5,
                                        signalValue=7,
                                        pvalue=8,
                                        qvalue=9)) #See readBroadPeak function
library(GenomicRanges)
# Lets take only chrosome 21
ctcf.peaks = ctcf.peaks[seqnames(ctcf.peaks) == 'chr21'] 
ctcf.peaks = ctcf.peaks[order(-ctcf.peaks$signalValue)]
summary(width(ctcf.peaks))
# Some peaks are extremly long, so lets trim them to 1kb around center
ctcf.peaks = resize(ctcf.peaks, width=1000, fix='center')

# Read transcripts from the refseq bed file 
# within 4 kp distance from the longest transcript of each gene
refseq.path <- "./refseq.hg19.bed" # I downloaded this file from the UCSC website
transcriptFeat=readTranscriptFeatures(refseq.path, up.flank = 2000, down.flank = 2000) # it's a GRangesList object
# Statistics of overlap of our CTCF peaks with gene structures like introns, exons etc.
ann.ctcf <- annotateWithGeneParts(ctcf.peaks, transcriptFeat ) 
# Calculate percentage of CTCF peaks overlapping with annotation
getTargetAnnotationStats(ann.ctcf, percentage=TRUE, precedence=TRUE)
# Plot it as pie chart
plotTargetAnnotation(ann.ctcf, main="CTCF in H1hesc", cex.legend = .7)

# Let's check overlap of ctcf peaks with annotation from chromHMM 
chrHMM.url="http://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeBroadHmm/wgEncodeBroadHmmH1hescHMM.bed.gz"
chrHMM=readBed(chrHMM.url)
chrHMM.list=GenomicRanges::split(chrHMM, chrHMM$name, drop=TRUE)
# Read more windows
broadpeak.files = list.files(genomationDataPath, full.names=TRUE, pattern='.broadPeak.gz$')
broadpeak.list = GRangesList( lapply(broadpeak.files, readGeneric) )
names(broadpeak.list) <- c("Ctcf","P300","Suz12", "Rad21")
peak2ann.l=annotateWithFeatures(broadpeak.list, chrHMM.list)
# Percentage of target elements (Ctcf, P300, etc) overlapping with features (chrMM annotaion)
heatTargetAnnotation(peak2ann.l)

# Create matrix of coverage of reads from given bam files on Ctcf peaks
sm = ScoreMatrix(bam.files[1], ctcf.peaks, type='bam')
# Bin score matrix into 50 bins and take average within bins (besides mean you can set max or min)
smb = ScoreMatrixBin(bam.files[1], ctcf.peaks, bin.num=50, type='bam')
# Calculate many score matrices at once in parallel
sml = ScoreMatrixList(bam.files, ctcf.peaks, bin.num=50, type='bam', cores=2)
# Descriptions of file that contain info. about transcription factors
sampleInfo = read.table(system.file('extdata/SamplesInfo.txt',
                                    package='genomationData'),header=TRUE, sep='\t')
names(sml) = sampleInfo$sampleName[match(names(sml),sampleInfo$fileName)]
# Visualite score matrix
multiHeatMatrix(sml, xcoords=c(-500, 500))
heatMeta(sml, xcoords=c(-500, 500))

# Scale each ScoreMatrix in the ScoreMatrixList object, by rows and/or columns
# Due to large signal scale of rows of each element in the ScoreMatrixList lets scale them
sml.scaled = scaleScoreMatrixList(sml)
multiHeatMatrix(sml.scaled, xcoords=c(-500, 500))
heatMeta(sml.scaled, xcoords=c(-500, 500))

# k-means with k=2
cl1 <- function(x) kmeans(x, centers=2)$cluster
multiHeatMatrix(sml.scaled, xcoords=c(-500, 500), clustfun = cl1)

# hierarchical clustering with Ward's method for agglomeration into 2 clusters
cl2 <- function(x) cutree(hclust(dist(x), method="ward"), k=2)
multiHeatMatrix(sml.scaled, xcoords=c(-500, 500), clustfun = cl2)

# Defining which matrices are used for clustering: “clust.matrix” in multiHeatMatrix
multiHeatMatrix(sml.scaled, xcoords=c(-500, 500), clustfun = cl1, clust.matrix = 1)

# Central tendencies in line plots: centralTend in plotMeta
# CentralTend arg can be median or mean
plotMeta(mat=sml.scaled, profile.names=names(sml.scaled),
         xcoords=c(-500, 500),
         winsorize=c(0,99),
         centralTend="mean")

plotMeta(mat=sml.scaled, profile.names=names(sml.scaled),
         xcoords=c(-500, 500),
         winsorize=c(0,99),
         centralTend="mean",  
         smoothfun=function(x) stats::smooth.spline(x, spar=0.5))

# Dispersion shows dispersion interval bands around central tendency
plotMeta(mat=sml, profile.names=names(sml),
         xcoords=c(-500, 500),
         winsorize=c(0,99),
         centralTend="mean",
         smoothfun=function(x) stats::smooth.spline(x, spar=0.5),
         dispersion="se", lwd=4)

# Calculate scores that correspond to k-mer or PWM matrix occurence: patternMatrix function
# Ctcf motif from the JASPAR database
ctcf.pfm = matrix(as.integer(c(87,167,281,56,8,744,40,107,851,5,333,54,12,56,104,372,82,117,402, 
                               291,145,49,800,903,13,528,433,11,0,3,12,0,8,733,13,482,322,181, 
                               76,414,449,21,0,65,334,48,32,903,566,504,890,775,5,507,307,73,266, 
                               459,187,134,36,2,91,11,324,18,3,9,341,8,71,67,17,37,396,59)), 
                  ncol=19,byrow=TRUE)
rownames(ctcf.pfm) <- c("A","C","G","T")

prior.params = c(A=0.25, C=0.25, G=0.25, T=0.25)
priorProbs = prior.params/sum(prior.params)
postProbs = t( t(ctcf.pfm + prior.params)/(colSums(ctcf.pfm)+sum(prior.params)) )
ctcf.pwm = Biostrings::unitScale(log2(postProbs/priorProbs))

library(BSgenome.Hsapiens.UCSC.hg19)
hg19 = BSgenome.Hsapiens.UCSC.hg19

p = patternMatrix(pattern=ctcf.pwm, windows=ctcf.peaks, genome=hg19, min.score=0.8)

# Visualization of the patternMatrix (here as ScoreMatrix object) can be 
# done by using i.e. heatMatrix, heatMeta or plotMeta functions.
heatMatrix(p, xcoords=c(-500, 500), main="Ctcf motif")
plotMeta(mat=p, xcoords=c(-500, 500), smoothfun=function(x) stats::lowess(x, f = 1/10), 
         line.col="red", main="Ctcf motif")



