\name{summary.cghFLasso}
\alias{summary.cghFLasso}
\alias{print.summary.cghFLasso}
\title{Summarizing gain/loss calls by cghFLasso}
\description{ 'summary' method for class \code{cghFLasso}.}
   
\usage{
         summary.cghFLasso(object, index,...)
         print.summary.cghFLasso(x,...)
}

\arguments{
    \item{object}{an object of class \code{cghFLasso} (returned by function \code{cghFLasso}).}
    \item{index}{numeric vector specifying which arrays to plot.} 
    \item{x}{an object of class \code{summary.cghFLasso}, usually, a result of a call to \code{summary.cghFLasso}.}
    \item{...}{further arguments passed to or from other methods.}
}
        
\details{
     \code{summary.cghFLasso} summarize the gain/loss calls for a group of CGH arrays. It reports the consensus counts of alterations
      as well as the corresponding FDRs for each gene/clones. It also returns the sample percentage of amplificaiton and deletion at each gene/cloens.
 
     \code{print.summary.cghFLasso} outputs a matrix summarizing the gain/loss calls across all the samples 
      for each chromosome. It returns the maximum sample percentage of amplificaiton and deletion on each chromosome respectively.
      It also reports the genome percentage of alteration for each chromosome (proportion of genes having consensus FDR smaller than 0.05 on the chromosome).  
}

\value{
 ans<-list(ConsensusCount=ConsensusCount, CC.FDR=CC.FDR, 
           Amp.ConCount=Amp.CC, Del.ConCount=Del.CC, 
           chrom.summary=chrom.summary, sample.summary=sample.summary,
           Esti.copy=Esti.copy, chromosome=chromosome, nucposi=nucposi)
    \item{ConsensusCount}{numeric vector with the same length as each CGH array. It reports the number of samples showing copy number alteration at each gene/clones}
    \item{CC.FDR}{numeric vector of the same length of \code{ConsensusCount}. It is the estiamted probability of observing 
                   the same or higher consensus count by random chance.}
    \item{Amp.ConCount}{numeric vector of the same length of \code{ConsensusCount}. 
                         It reports the number of samples showing amplificaiton at each gene/clones.}
    \item{Del.ConCount}{numeric vector of the same length of \code{ConsensusCount}. 
                         It reports the number of samples showing deletion at each gene/clones.}
    \item{chrom.summary}{matrix with four columns. The first column is the chromosome number. 
                         The second and third columns are the maximum sample percentages of amplificaitons and deletions 
                          on each chromosome respectively.
                         The forth column represents the genome percentage of alterations on each chromosome.}
    \item{sample.summary}{numeric vector with the same lenght of sample numbers. It reports the overall alteration percentage of each sample.}
    \item{Esti.copy}{numeric matrix showing estiamted copy numbers of all genes/clones of all samples.}
    \item{chromosome}{numeric vector with the same length as each CGH array. It's the chromosome number of each gene/clone.}
    \item{nucposi}{numeric vector with the same length as each CGH array. It's the nucletide position of each gene/clone.}
  }

\references{
R. Tibshirani, M. Saunders, S. Rosset, J. Zhu and K. Knight (2004) `Sparsity and smoothness via the fused lasso', 
J. Royal. Statist. Soc. B. (In press), available at http://www-stat.stanford.edu/~tibs/research.html.

P. Wang, Y. Kim, J. Pollack, B. Narasimhan and R. Tibshirani (2005)
`A method for calling gains and losses in array CGH data', 
Biostatistics 2005, 6: 45-58, available at http://www-stat.stanford.edu/~wp57/CGH-Miner/

R. Tibshirani and P. Wang (2007) `Spatial smoothing and hot spot detection using the Fused Lasso',
Biostatistics (In press), available at http://www-stat.stanford.edu/~tibs/research.html.

J. Friedman, T. Hastie. R. Tibshirani (2007) `Pathwise coordinate optimization and  the fused lasso'.
}

\author{N. A. Johnson, R. Tibshirani and P. Wang}

\keyword{methods}

\examples{



library(cghFLasso)
data(CGH)

#############
### Example 1: Process one chromosome vector without using normal references.

CGH.FL.obj1<-cghFLasso(CGH$GBM.y)
plot(CGH.FL.obj1, index=1, type="Lines")

#############
### Example 2: Process a group of CGH arrays and use normal reference arrays.

Normal.FL<-cghFLasso.ref(CGH$NormalArray,  chromosome=CGH$chromosome)
Disease.FL<-cghFLasso(CGH$DiseaseArray, chromosome=CGH$chromosome, nucleotide.position=CGH$nucposition, FL.norm=Normal.FL, FDR=0.01)

###  Plot for the first arrays
i<-1
plot(Disease.FL, index=i, type="Single")
title(main=paste("Plot for the ", i ,"th BAC array", sep=""))

### Consensus plot
plot(Disease.FL, index=1:4, type="Consensus")
title(main="Consensus Plot for 4 BAC arrays")

### Plot all arrays
plot(Disease.FL, index=1:4, type="All")
title(main="Plot for all 4 arrays")

### Report and output
report<-summary(Disease.FL, index=1:4)
print(report)
output.cghFLasso(report, file="CGH.FL.output.txt")
}