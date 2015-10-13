# 'FLS' R group - October 2015
# Patterns and features on genomes in R: Biostrings, BSgenome and GenomicRanges.
# Dave Gerrard
# Biostrings, BSgenomes and GenomicRanges.

# PRE-REQUISITES --------------------------------
# Multiple bioconductor packages but should get them all by doing this. 
## try http if https is not available
#source("https://bioconductor.org/biocLite.R")
#biocLite("BSgenome.Scerevisiae.UCSC.sacCer3")     # 'should' also install Biostrings, BSgenomes and GenomicRanges if not already installed.

# setwd("/folder/with/data/and/script")

# The Plan. --------------------------------------
# Find instances of a DNA-sequence motif in a genome that also have evidence of binding from a ChIP-seq experiment.
# 1. Biostrings
# 2. BSgenome
# 3. GenomicRanges

# Introduce Biostrings objects, reading (DNA) sequence, searching for patterns   ----------------------

library(Biostrings)
# http://bioconductor.org/packages/release/bioc/vignettes/Biostrings/inst/doc/Biostrings2Classes.pdf
lsf.str("package:Biostrings")   # recent tip on @RLangTip

b <- BString("I am a BString object")
b 
length(b)
is(b)                                        # Xstrings
d <- DNAString("I am a BString object")      # should error
d <- DNAString("TTGAAAA-CTC-N")
is(d)
# what letter allowed?   (plus gaps)
IUPAC_CODE_MAP

d[1:10]
subseq(d, 3,8)                               # better for large strings.
reverseComplement(d)

# load example DNA data from Biostrings package
filepath1 <- system.file("extdata", "someORF.fa", package="Biostrings")  # go and have a look, this is just a fasta file.
filepath1
fasta.info(filepath1, seqtype="DNA")
x1 <- readDNAStringSet(filepath1)
x1
x1[1:2]
subseq(x1[1:2], start=c(1,2), end=10)  # grab some fragments.

matchPattern('TTGTAAATATATCTT', x1)
matchPattern('TTGTAAATATATCTT', x1[[1]])      # there are not vmatchPattern functions for every occasion (yet)

matches <- vmatchPattern('TTGTAAATATATCTT', x1)
matches
matches[[1]]
# useful to get counts
counts <- vcountPattern('TTGTAAATATATCTT', x1)


# Let's use a real motif
# Ste12 motif http://www.nature.com/nature/journal/v464/n7292/fig_tab/nature08934_F3.html
# ideally get the pwm and use matchPWM
# Here, we'll just use a single sequence representation
pattern.core <- "GAAAC"
pattern.long <- "TGAAACR"  # R = puRine (A,G)
IUPAC_CODE_MAP

vmatchPattern(pattern.core, x1)
vmatchPattern(pattern.long, x1)         # no hits?
vcountPattern(pattern.long, x1)
vcountPattern(pattern.long, x1, fixed=F)  # to use degenerate base codes
vmatchPattern(pattern.long, x1, fixed=F)

# Introduce BSgenome packages. load sacCer3 as BSgenome  -----------------
library(BSgenome.Scerevisiae.UCSC.sacCer3)    # note the other packages being loaded.
available.genomes()
# ask me another time how to make your own BSgenome

# search for motif in a genome.
genome <- BSgenome.Scerevisiae.UCSC.sacCer3   # for convenience
genome
# accessor functions:-
seqnames(genome)
seqlengths(genome)
barplot(seqlengths(genome), horiz=T, las=1)   # today's only plot?
summary(genome)
organism(genome)
provider(genome)    # ... etc

# searching one chromosome
thisChrom <- genome[["chrI"]]  # when is it loaded into memory?
thisChrom    # with some genomes (e.g. Human) you will just see NNNNNNNNNNNNNNNNNNNN...NNNNNNNNNNNNNNNNNNNN - unsequenced telomeres.
#matchPattern("TTCCCTTC", thisChrom, fixed=F)
matchPattern(pattern.long, thisChrom, fixed=F)
matchPattern(pattern.core, thisChrom, fixed=F)      # Views objects...

# To match the whole genome use vmatchPattern.  (but no vmatchPWM, so may need to iterate through chroms).
tf.hits <- vmatchPattern(pattern.long, genome, fixed=F)  # 
tf.hits 
is(tf.hits )   # GenomicRanges

# load peak regions for  chip-seq data that should match motif.  convert to gRanges.

# e.g. Ste12 in yeast  http://www.nature.com/nature/journal/v464/n7292/full/nature08934.html
# files from http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE19635 
# strain is S288c  need SacCer3  
# Introduce 'hits' objects. Load data and create GRanges object
peaks.raw <- read.delim("GSE19635_s96a_peaks.txt", comment.char = "#")  # edit the location to work on your machine.
head(peaks.raw)
dim(peaks.raw)
# chromosome field does not match exactly with seqnames(genome)
levels(peaks.raw$chr)                             # no chrM
# todo
# quick and risky
old.counts <- table(peaks.raw$chr)               # store this for simple check after conversion
levels(peaks.raw$chr) <- paste("chr", as.roman(1:16), sep="")
head(peaks.raw)
# might be better to individually remake each level. 
old.counts 
table(peaks.raw$chr)


# make GenomicRanges object.
head(peaks.raw)
peaks.GR <- GRanges(peaks.raw$chr, IRanges(peaks.raw$start, peaks.raw$end))   # use of IRanges() within GRanges()
peaks.GR

# findOverlaps and other 'set' operations in GenomicRanges.
# as we are interested in finding the motifs that are bound, we'll make the tf.hits the query.
tf.hits
ov <- findOverlaps(tf.hits, peaks.GR)
ov                      # a pair of indices
ova <- overlapsAny(tf.hits, peaks.GR)
head(ova, n=100)

tf.hits[ov@queryHits[1]]
peaks.GR[ov@subjectHits[1]]
width(peaks.GR[ov@subjectHits])

this.peak <- peaks.GR[ov@subjectHits[1]]
this.peak.seq <- getSeq(genome, names=seqnames(this.peak), start=start(this.peak), end=end(this.peak))

# and just because I can.
pairwiseAlignment(pattern = pattern.long, subject = this.peak.seq, type="local-global")

# how many ChIP-peaks contain at least one motif
sum(overlapsAny(peaks.GR,tf.hits))
length(peaks.GR)

# get the sequences from those peaks
motifPeaks <- peaks.GR[unique(ov@subjectHits)]    # there may be duplicates.  inspect 'ov'
motifPeak.seqs <- getSeq(genome, names=seqnames(motifPeaks), start=start(motifPeaks), end=end(motifPeaks))
names(motifPeak.seqs)        # missing useful names
names(motifPeak.seqs) <- paste0(seqnames(motifPeaks), ":", start(motifPeaks), "-", end(motifPeaks))
writeXStringSet(motifPeak.seqs, file="motifPeaks.fasta")


# from here....
# rtracklayer
# motif refinement/discovery
# motif co-occurrence
# permutation testing
# gene associations
# gene-set enrichment
# SNP assocationss
# export data  head(as.data.frame(peaks.GR))

# tip
?reduce  # this and other clever set operations.
# not covered
# masks, matchPDict(), letterfrequency()

# if time
# packrat::status()  # this probably won't work on your machine..
