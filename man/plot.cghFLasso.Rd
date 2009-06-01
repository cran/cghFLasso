\name{plot.cghFLasso}
\alias{plot.cghFLasso}
\title{A function to plot gain/loss calls on CGH arrays}
\description{A function to plot gain/loss calls on CGH arrays}
   
\usage{
plot.cghFLasso(x, index, type="All", centro=NULL, ...)
}
\arguments{
    \item{x}{an object of class \code{cghFLasso} (returned by function \code{cghFLasso}).}
    \item{index}{numeric vector specifying which arrays to plot.} 
    \item{type}{a character specifying the plot type. It should be one of "Lines", "Single", "Consensus", or "All". See details for more information.}
    \item{centro}{numeric vector specifying the centromere positions. If missing, the default centromere value of human genome will be used.}
    \item{...}{further arguments passed to or from other methods.}
}
        
\details{
    \code{plot.cghFLasso} provides several summary plots for \code{cghFLasso} result object. If \code{type="Lines"}, 
     it plots the raw CGH measurements and the gain/loss results for the one selected array along the genome (no seperation of chromosomes). If \code{type="Single"}, 
     it plots the raw CGH measurements and the gain/loss results for the one selected array along the genome (different chromosomes are ploted in 
     different lines). If \code{type="Consensus"}, 
     it plots the percentages of the gain/loss calls across a group of samples along the genome. If \code{type="All"}, it plots the results of 
     all selected arrays in one figure.      
}

\value{
     It returns the \code{type} of the plot.    
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