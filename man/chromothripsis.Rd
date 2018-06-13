\name{chromothripsis} 
\alias{chromothripsis} 
\title{CNomplexity::chromothripsis} 
\usage{chromothripsis(bedpeFile, size=3e7,
			gap=1e6,
			bedpeChromCol1=1,bedpePosCol1=2,direction1col=9,
			bedpeChromCol2=4,bedpePosCol2=5,direction2col=10,
			bedpeSampleCol=1,svclassCol=NULL,
			doParallel=FALSE,
			nCores = NULL,samplesToRun = NULL,
			chromsToRun = NULL,sepbedpe=TRUE,
			bedpeHead=FALSE,bedpeEnding=".brass.annot.bedpe.gz",
			breaksLimit=30,pThresh=0.8,
			cytoFile=NULL,
			diag=FALSE)}
\arguments{
\item{bedpeFile}{Directory of separate bedpe files, or filename of single bedpe file.}
\item{size}{Window size for sliding window across chomosomes.}
\item{gap}{Gap between consecutive sliding windows.}
\item{bedpeChromCol1}{Column of bedpe file with chromosome name for first fusion partner.}
\item{bedpePosCol1}{Column of bedpe file with position for first fusion partner.}
\item{direction1col}{Column of bedpe file with orientation for first fusion partner.}
\item{bedpeChromCol2}{Column of bedpe file with chromosome name for second fusion partner.}
\item{bedpePosCol2}{Column of bedpe file with position for second fusion partner.}
\item{direction2col}{Column of bedpe file with orientation for second fusion partner.}
\item{bedpeSampleCol}{Column of bedpe file with sample name.}
\item{svclassCol}{Column of structural variation class e.g. deletion, inversion, tandem-duplication or translocation}
\item{doParallel}{Whether to run in parallel}
\item{nCores}{Number of cores if \code{doParallel=TRUE}.}
\item{samplesToRun}{Character vector with sample names to run analysis for.}
\item{chromsToRun}{Character vector with chromosome names to run analysis for.}
\item{sepbedpe}{Whether bedpe files are separate or concatenated.}
\item{bedpeHead}{Whether bedpe file has a header.}
\item{bedpeEnding}{File ending for separate bedpe files.}
\item{breaksLimit}{Minimum number of breakpoints on a chromsome to consider chromothripsis. Defaults to 30.}
\item{pThresh}{P value threshold for test of randomness (random order of breakpoints and random joins). Defaults to 0.8.}
\item{cytoFile}{File with starts and ends of chromosome bands. If \code{NULL} defaults to a file supplied with CNomplexity.}
\item{diag}{Boolean. Whether to return diagnostic outputs.}
}
\value{
List of samples and chromosmores, with each element being a matrix with Sample name, start positon of chromothriptic region and end position of chromothriptic region. NULL is no chromothriptic regions in the chromosome.
} 
\description{
Identify chromothriptic regions. 
}
\details{
  Implements three tests from Korbel & Campbell (2013).
}

