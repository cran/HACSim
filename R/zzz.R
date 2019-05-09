envr <- NULL 

.onLoad <- function(...) {
  envr <<- new.env()  # when package is loaded, create new environment to store needed variables 
}

.onAttach <- function(...) {
  packageStartupMessage("This is HACSim 1.0.0 \n
Type ?HACHypothetical to see how to set up objects to run \n simulations of haplotype accumulation for hypothetical species \n 
Type ?HACReal to see how to set up objects to run \n simulations of haplotype accumulation for real species \n 
Type ?HAC.sim to see how run simulations of haplotype \n accumulation curves")
}