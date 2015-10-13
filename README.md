Email from Dave:

> For those who would like to follow along on their laptop in tomorrow's session, you will need to install some Bioconductor packages and bring the attached script file [BSgenome.BioStrings.bioC.R] and data file [GSE19635_s96a_peaks.txt].
>
> To install the Bioconductor packages, connect to the internet, open your R prompt or Rstudio and type (copy/paste) these two lines:-
>

> ```
source("https://bioconductor.org/biocLite.R")
biocLite("BSgenome.Scerevisiae.UCSC.sacCer3")    
```

> This will install a yeast genome as a package and should also install Biostrings, BSgenomes and GenomicRanges if not already installed. If you've not used Bioconductor before, this may take a few minutes to download and install everything.

The 2 attached files are uploaded in this repository.

The previous email:

> Hi everyone,
> 
> We will be meeting tomorrow (Tue 13th) in the Michael Smith Building, room A.3025, 11AM to 12PM.
> 
> Dave Gerrard will be presenting:
> Patterns and features on genomes in R: Biostrings, BSgenome and GenomicRanges.
> 
> "I will give a brief introduction of using some key genomics packages from Bioconductor.  I'll give some simple examples of reading in sequences, searching for patterns (motifs), getting and using published genomes and finding co-occurrence of features (e.g. motifs and ChIP-seq peaks) on the genome."
> 
> Bring a mug for tea and coffee.
> 
> The timetable for the next few weeks can be found [https://github.com/flsrgroup/Home](https://github.com/flsrgroup/Home)
