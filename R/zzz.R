### Functions needed when loading HACSim ###

envr <- NULL

.onLoad <- function(...) {
  envr <<- new.env() # when package is loaded, create new environment to store needed variables
}

.onAttach <- function(...) {
  packageStartupMessage("This is HACSim 1.0.6 \n
Type ?HACHypothetical to see how to set up objects to run \n simulations of haplotype accumulation for hypothetical species \n 
Type ?HACReal to see how to set up objects to run \n simulations of haplotype accumulation for real species \n 
Type ?HAC.simrep to see how to run simulations of haplotype \n accumulation curves \n
Type ?sim.seqs to see how to simulate DNA sequences according to \n various models of nucleotide substitution \n
Type ?launchApp to see how to run the HACSim Shiny app")
}
