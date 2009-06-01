\name{CGH}
\alias{CGH}
\title{Example CGH Array Data}
\description{A list of four components (an example CGH array data for package \code{cghFLasso})}

\format{
Both \code{NormalArray} and \code{DiseaseArray} are muneric matrix of 2270 rows and 3 columns: each row corresponds to one gene/clone, 
and each column corresponds to one CGH array. The value of each entry is the log fluorescence
ratio resulted from the CGH experiment. The order of the genes/clones in the rows is the same as the order of the
genes/clones on the genome. \code{chromosome} and \code{nucposition} provide chromosome number and nucleotide 
position for each gene/clone. \\
\code{GBM.y} is a numeric vertor of length 990, providing CGH measurement of one pseudo chromosome.
}

\references{
P. Wang, Y. Kim, J. Pollack, B. Narasimhan and R. Tibshirani, 
"A method for calling gains and losses in array CGH data", Biostatistics 2005, 6: 45-58, 
available at http://www-stat.stanford.edu/~wp57/CGH-Miner/

R. Tibshirani and P. Wang (2007) `Spatial smoothing and hot spot detection using the Fused Lasso',
Biostatistics (In press), available at http://www-stat.stanford.edu/~tibs/research.html.
}

\usage{data(CGH)}

\keyword{datasets}

