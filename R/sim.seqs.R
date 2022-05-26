### Simulate DNA sequences according to various models of nucleotide evolution ###
# Experimental functionality at this point 

#####

sim.seqs <- function(num.seqs,
                     num.haps, 
                     length.seqs,
                     count.haps,
                     nucl.freqs,
                     codon.tbl = c("standard", "vertebrate mitochondrial", "invertebrate mitochondrial"),
                     subst.model = c("JC69", "K80", "F81", "HKY85"),
                     mu.rate,
                     transi.rate,
                     transv.rate) {
  
  # Error messages 
  
  if (num.seqs < num.haps) {
    stop("Number of sequences must be greater than or equal to number of haplotypes")
  }
  
  if (num.seqs == 1) {
    stop("Number of sequences must be greater than 1")
  }
  
  if (num.haps == 1) {
    stop("Number of haplotypes must be greater than 1")
  }
  
  if (!isTRUE(all.equal(1, sum(nucl.freqs), tolerance = .Machine$double.eps^0.25))) {
    stop("Nucleotide frequencies must sum to 1")
  }
  
  if (length(count.haps) != num.haps) {
    stop("count.haps must have num.haps elements")
  }
  
  if (sum(count.haps) != num.seqs) {
    stop("Sum of haplotype counts must equal the number of sequences")
  }
  
  if (((subst.model == "JC69") || (subst.model == "K80")) && (all(nucl.freqs != 0.25)))  {
    stop("Both JC69 and K80 DNA substitution models require all nucleotide frequencies to equal 0.25")
  }
  
  if (((subst.model == "F81") || (subst.model == "HKY85")) && (all(nucl.freqs == 0.25)))  {
    stop("Both F81 and HKY85 DNA substitution models require all nucleotide frequencies to differ from 0.25")
  }
  
  if (((subst.model == "JC69") || (subst.model == "F81")) && (is.null(mu.rate) == TRUE))  {
    stop("Both JC69 and F81 DNA substitution models require the overall nucleotide mutation rate to be specified")
  }
  
  if (((subst.model == "K80") || (subst.model == "HKY85")) && ((is.null(transi.rate) == TRUE) || ((is.null(transv.rate) == TRUE))))  {
    stop("Both K80 and HKY85 DNA substitution models require each of the nucleotide transition rate and nucleotide transversion rate to be specified")
  }
  
  # DNA alphabet
  
  nucl <- as.DNAbin(c("a", "c", "g", "t")) # A, C, G, T
  
  # Generate a single random DNA sequence of given length according to nucleotide frequency distribution
  # This is the "seed" sequence from which all other sequences are generated
  seqs <- sample(nucl, size = length.seqs, replace = TRUE, prob = nucl.freqs)
  
  # DNA codons
  
  pos1 <- c("a", "c", "g", "t")
  pos2 <- c("a", "c", "g", "t")
  pos3 <- c("a", "c", "g", "t")
  
  codons <- expand.grid(pos1, pos2, pos3)
  codons <- paste0(codons$Var1, codons$Var2, codons$Var3) # 64 possible codons
  
  # Exclude stop codons according to desired genetic code
  
  if (codon.tbl == "standard") {
    stop.codons <- c("taa", "tag", "tga")
  } else if (codon.tbl == "vertebrate mitochondrial") {
    stop.codons <- c("aga", "agg", "taa", "tag")
  } else {
    # invertebrate mitochondrial
    stop.codons <- c("taa", "tag")
  }
  
  # Allowable codons
  codons <- codons[!codons %in% stop.codons]
  
  regx <- paste0("(", paste(stop.codons, collapse = ")|("), ")")
  s <- paste(seqs, collapse = "")
  seqs <- gsub(regx, paste(sample(nucl, size = 3, replace = TRUE, prob = nucl.freqs), collapse = ""), s)
  seqs <- strsplit(seqs, "")
  
  seqs <- as.matrix(as.DNAbin(seqs))
  
  # DNA substitution models

  if ((subst.model == "JC69") || (subst.model == "F81")) {
    
    # Define mutations - all mutations are equally likely 
    
    mu.set <- list('a' = as.DNAbin('c'),
                   'a' = as.DNAbin('g'),
                   'a' = as.DNAbin('t'),
                   'c' = as.DNAbin('a'),
                   'c' = as.DNAbin('g'),
                   'c' = as.DNAbin('t'),
                   'g' = as.DNAbin('a'),
                   'g' = as.DNAbin('c'),
                   'g' = as.DNAbin('t'),
                   't' = as.DNAbin('a'),
                   't' = as.DNAbin('c'),
                   't' = as.DNAbin('g'))
    
    muts <- function(seqs) {
      unlist(mu.set[as.character(seqs)])
    }
    
  } else {
    # Define transitions and transversion mutations
    
    transi.set <- list('a' = as.DNAbin('g'), 
                       'c' = as.DNAbin('t'),
                       'g' = as.DNAbin('a'), 
                       't' = as.DNAbin('c'))
    transv.set <- list('a' = as.DNAbin(c('c', 't')),
                       'c' = as.DNAbin(c('a', 'g')),
                       'g' = as.DNAbin(c('c', 't')), 
                       't' = as.DNAbin(c('a', 'g')))
    
    transi <- function(seqs) {
      unlist(transi.set[as.character(seqs)])
    }
    
    transv <- function(seqs) {
      # randomly sample a transversion mutation
      sapply(transv.set[as.character(seqs)], FUN = sample, 1)
    }
  }
  
  duplicate.seq <- function(seqs) {
    if ((subst.model == "JC69") || (subst.model == "F81")) {
      # calculate number of sequence substitutions given desired mutation rate according to binomial distribution
      num.muts <- rbinom(n = 1, size = length.seqs, prob = mu.rate) 
      if (num.muts > 0) {
        # randomly place mutations along the DNA sequence
        idx <- sample(length.seqs, size = num.muts, replace = FALSE) # segregating sites
        seqs[idx] <- muts(seqs[idx])
      } 
    } else {
      # calculate number of sequence transitions given desired mutation rate according to binomial distribution
      num.transi <- rbinom(n = 1, size = length.seqs, prob = transi.rate) # total number of transitions
      if (num.transi > 0) {
        # randomly place transitions mutations along the DNA sequence
        idx <- sample(length.seqs, size = num.transi, replace = FALSE) # segregating sites
        seqs[idx] <- transi(seqs[idx])
      }
      
      # calculate number of sequence transitions given desired mutation rate according to binomial distribution
      num.transv <- rbinom(n = 1, size = length.seqs, prob = transv.rate) # total number of transversions
      if (num.transv > 0) {
        # randomly place transitions mutations along the DNA sequence
        idx <- sample(length.seqs, size = num.transv, replace = FALSE) # segregating sites
        seqs[idx] <- transv(seqs[idx])
      }
    }
    seqs
  }
  
  # Generate num.seqs random DNA sequences from "seed"sequence based on desired nucleotide substitution model
  seqs <- matrix(replicate(num.seqs, duplicate.seq(seqs)), byrow = TRUE, nrow = num.seqs)
  
  # Randomly sample generated DNA sequences based on haplotype frequency distribution
  seqs <- seqs[sample(num.haps, size = num.seqs, replace = TRUE, prob = count.haps), ]

  # convert to DNAbin object
  class(seqs) <- "DNAbin" 
  
  # write DNA sequences to a FASTA file called 'simseqs.fas'
  write.dna(seqs, file = "simseqs.fas", format = "fasta")
  
} 
