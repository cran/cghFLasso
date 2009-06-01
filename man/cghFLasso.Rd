\name{cghFLasso}
\alias{cghFLasso}
\title{A function to call alteration for CGH arrays using the fused lasso regression}
\description{A function to call alteration for CGH arrays using the fused lasso regresssion}
   
\usage{
cghFLasso(CGH.Array, chromosome=NULL, nucleotide.position=NULL, FL.norm=NULL, FDR=NULL, filter=NULL, missing.PlugIn=TRUE, smooth.size=5)
}

\arguments{
    \item{CGH.Array}{numeric vector or matrix. It's the result of one or mutiiple CGH experiments.
                     Each column is the log2 ratios returned from one array experiment and
                     is ordered according to the gene/clones' position on the genome. Missing value should be coded as NA. 
                     }
    \item{chromosome}{numeric vector. Length should be the same as the row number of \code{CGH.Array}. 
                      It's the chromosome number of each gene/clone. If no value is specified, the arrays will be treated 
                      as one chromosome.} 
    \item{nucleotide.position}{numeric vector. Length should be the same as the row number of \code{CGH.Array}. 
                               It's the nucleotide position of each gene/clone. This information is used mainly for plot.
                               If no value is specified, the program will make genes/clones equally spaced on the genome.}  
    \item{FL.norm}{numeric vector or matrix. Smoothed result of the reference arrays. Set to NULL (default) if the reference arrays are not available. }
    \item{FDR}{numeric value (between 0 and 1). User can use this option to control False Discovery Rate of the results. 
               If not specified, the function will return the raw output of fused lasso regression on the target array.}
    \item{filter}{numeric vector. Length should be the same as the row number of \code{CGH.Array}. Each element takes value 1 or 0. Value 1 means that the 
                  corresponding gene/clone will be filter away in ahead of the analysis. If no \code{filter} is specified, the program will process all the 
                  genes/clones in the input arrays.}
    \item{missing.PlugIn}{Bollen value. If its value is TRUE, the missing values will be replaced with local averages.}
    \item{smooth.size}{numeric value, specifying the window size for carrying out the local average. The average is used only to 
                       replace the missing values.}
}
        
\details{
    \code{cghFLasso} calls copy number alterations for CGH arrays using the fused lasso regression.
    The dynamic programming algorithm developed by N.A.Johnson is used for the model fitting.  }

\value{
    an object of class \code{cghFLasso} with components
    \item{Esti.CopyN}{data matrix reporting the estimated DNA copy numbers for seleted genes/clones of all the samples.}
    \item{CGH.Array}{data matrix reporting the raw CGH array measurements for selected genes/clones of all the samples.}
    \item{chromosome, nucleotide.position}{numeric vectors reporting the chromosome numbers and the nucleotide positions of selected genes/clones. 
                      If \code{filter=NULL}, these are the same as the input \code{chromosome} and \code{nucleotide.position}.}
    \item{FDR}{fdr value specified by the user.}
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