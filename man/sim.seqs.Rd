\name{sim.seqs}
\alias{sim.seqs}

\title{Simulate DNA sequences according to DNA substitution models}

\description{Simulates DNA sequences according to various DNA substitution 
models:

\itemize{
  \item Jukes-Cantor (1969)
  \item Kimura (1980)
  \item Felsenstein (1981)
  \item Hasegawa-Kishino-Yano (1985)
}

}

\usage{sim.seqs(num.seqs, num.haps, length.seqs, count.haps, nucl.freqs, 
codon.tbl = c("standard", "vertebrate mitochondrial", 
"invertebrate mitochondrial"), subst.model = c("JC69", "K80", "F81", "HKY85"), 
mu.rate, transi.rate, transv.rate)}

\arguments{\item{num.seqs}{Number of simulated DNA sequences}
  \item{num.haps}{Number of simulated unique species' haplotypes}
    \item{length.seqs}{Basepair length of DNA sequences}
    \item{count.haps}{Haplotype frequency distribution vector}
    \item{nucl.freqs}{Nucleotide frequency distribution vector of A, C, G, and 
    T respectively}
    \item{codon.tbl}{Codon table}
    \item{subst.model}{Model of DNA substitution}
    \item{mu.rate}{Overall nucleotide mutation rate/site/generation}
    \item{transi.rate}{Nucleotide transition rate/site/generation}
    \item{transv.rate}{Nucleotide transversion rate/site/generation}
  }
    
    \value{A FASTA file of DNA sequences}
    
\note{\code{num.seqs} must be greater than or equal to \code{num.haps}.

  Both \code{num.seqs} and \code{num.haps} must be greater than 1.

  \code{nucl.freqs} must have a length of four and its elements must sum to 1.
  
  \code{count.haps} must have a length of \code{num.haps} and its elements must 
  sum to \code{num.seqs}. 
  
  \code{subst.model} must be one of "JC69" (Jukes Cantor corrected p-distance), 
  "K80" (Kimura-2-Parameter (K2P), "F81" (Felenstein) or "HKY85"
  (Hasegawa-Kishino-Yano)
  
  \code{mu.rate} must be specified for both "JC69" and "F81" models
  
  \code{transi.rate} and \code{transv.rate} must be specified for both "K80" 
  and "HKY85" models
  
  All elements \code{nucl.freqs} must be equal to 0.25 when \code{subst.model} 
  is either "JC69" or "K80" 
  
  All elements \code{nucl.freqs} must differ from 0.25 when \code{subst.model} 
  is either "F81" or "HKY85"
  
  }

\examples{

\dontrun{

# Simulate DNA sequences from the 5'-COI DNA barcode region under a Jukes 
# Cantor nucleotide substitution model

num.seqs <- 100 # number of DNA sequences
num.haps <- 10 # number of haplotypes
length.seqs <- 658 # length of DNA sequences
count.haps <- c(60, rep(10, 2), rep(5, 2), rep(1, 5)) # haplotype frequency distribution
nucl.freqs <- rep(0.25, 4) # nucleotide frequency distribution
codon.tbl <- "vertebrate mitochondrial"
subst.model <- "JC69" # desired nucleotide substitution model
mu.rate <- 1e-3 # mutation rate
transi.rate <- NULL # transition rate
transv.rate <- NULL # transversion rate

sim.seqs(num.seqs = num.seqs, num.haps = num.haps, length.seqs = length.seqs, 
count.haps = count.haps, nucl.freqs = nucl.freqs, subst.model = subst.model, 
codon.tbl = codon.tbl, transi.rate = transi.rate, transv.rate = transv.rate)

}

}