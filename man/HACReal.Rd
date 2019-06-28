\name{HACReal}
\alias{HACReal}

\title{Set up an object to simulate haplotype accumulation curves for a real species}

\description{Helper function which creates an object containing necessary information to run a simulation of haplotype accumulation for a real species of interest}

\usage{HACReal(perms = 10000, p = 0.95, conf.level = 0.95, 
subsample = FALSE, prop = NULL, progress = TRUE, filename = NULL)}

\arguments{\item{perms}{Number of permutations (replications)}
\item{p}{Proportion of haplotypes to recover}
\item{conf.level}{Desired confidence level for graphical output and interval estimation}
\item{subsample}{Is a subsample of DNA sequences desired?}
\item{prop}{If subsample = TRUE, the proportion of DNA sequences to subsample}
\item{progress}{Should iteration output be printed to the R console?}
\item{filename}{Name of file where simulation results are to be saved}
}

\value{An object with 13 elements that can be passed to \code{HAC.simrep}
}

\examples{
    ## Simulate real species ##
    # outputs file called "output.csv"
    HACSObj <- HACReal(filename = "output") 
    
    ## Simulate real species - subsampling ##
    # subsamples 25\% of DNA sequences
    HACSObj <- HACReal(subsample = TRUE, prop = 0.25, conf.level = 0.95, 
    filename = "output") 
    
    ## Simulate real species and all parameters changed - subsampling ##
    HACSObj <- HACReal(perms = 10000, p = 0.90, subsample = TRUE, 
    prop = 0.15, conf.level = 0.99, filename = "output")
}
