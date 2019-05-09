\name{HAC.simrep}
\alias{HAC.simrep}

\title{Run a simulation of haplotype accumulation curves for hypothetical or real species}

\description{Runs the \code{HACSim} algorithm to iteratively extrapolate haplotype 
accumulation curves to determine likely specimen sample sizes for hypothetical or real species
}

\usage{HAC.simrep(HACSObject)}

\arguments{\item{HACSObject}{object containing desired simulation parameters}
}

\value{Iteration results are outputted to the console and graphs displayed in plot window. Additionally, iteration results are optionally saved to a CSV file. Subsampled DNA sequences are saved to a FASTA file.}

\note{When simulating real species via \code{HACReal(...)}, a pop-up window will appear prompting the user to select an intraspecific FASTA file of aligned/trimmed DNA sequences. The alignment must not contain missing or ambiguous nucleotides (i.e., it should only contain A, C, G or T); otherwise, haplotype diversity may be overestimated. Excluding sequences or alignment sites with missing/ambiguous data is an option.
}

\examples{

  ## Simulate hypothetical species ##
  
  N <- 100 # total number of sampled individuals
  Hstar <- 10 # total number of haplotypes
  probs <- rep(1/Hstar, Hstar) # equal haplotype frequency distribution
  
  HACSObj <- HACHypothetical(N = N, Hstar = Hstar , probs = probs, 
  filename = "output") # outputs a CSV file called "output.csv"
  
  ## Simulate hypothetical species - subsampling ##
  HACSObj <- HACHypothetical(N = N, Hstar = Hstar, probs = probs, 
  perms = 1000, p = 0.95, subsample = TRUE, prop = 0.25, 
  conf.level = 0.95, filename = "output")
  
  ## Simulate hypothetical species and all paramaters changed - subsampling ##
  HACSObj <- HACHypothetical(N = N, Hstar = Hstar, probs = probs, 
  perms = 10000, p = 0.90, subsample = TRUE, prop = 0.15, conf.level = 0.95, 
  filename = "output")
  
  HAC.simrep(HACSObj) # runs a simulation
  
  ## Simulate real species ##
  
  \donttest{
    ## Simulate real species ##
    # outputs file called "output.csv"
    HACSObj <- HACReal(filename = "output") 
    
    ## Simulate real species - subsampling ##
    HACSObj <- HACReal(subsample = TRUE, prop = 0.15, conf.level = 0.95, 
    filename = "output")
    
    ## Simulate real species and all parameters changed - subsampling ##
    HACSObj <- HACReal(perms = 10000, p = 0.90, subsample = TRUE, 
    prop = 0.15, conf.level = 0.99, filename = "output")
    
    # user prompted to select appropriate FASTA file
    HAC.simrep(HACSObj) 
    }
}