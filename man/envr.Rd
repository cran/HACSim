\name{envr}
\alias{envr}

\title{Simulation variable storage environment}

\description{\code{envr} is a new (initially empty) environment that is created 
when \code{HACSim} is loaded. 
}

\value{When a simulation is run via \code{HAC.simrep}, \code{envr} will contain 
26 elements as follows:

\item{ci.type}{Type of confidence interval to compute and plot. Defaul is 
\code{conf.type = "quantile"}.}
\item{conf.level}{The desired confidence level. Default is 
\code{conf.level = 0.95}.}
\item{d}{A dataframe with \code{Nstar - X} rows and five columns: specimens 
(specs), accumulated haplotypes (means), standard deviations (sds) and 
quantiles (both lower and upper)}
\item{df.out}{A dataframe with \code{iters} rows and six columns displaying 
"Measures of Sampling Closeness".}
\item{filename}{The name of the file where results are to be saved. Default 
is NULL.}
\item{Hstar}{Number of unique species' haplotypes}
\item{input.seqs}{Should DNA sequences be inputted? Default is FALSE.}
\item{iters}{The number of iterations required to reach convergence}
\item{N}{The starting sample size used to initialize the algorithm}
\item{Nstar}{The final (extrapolated) sample size}
\item{Nstar.high}{The upper endpoint of the desired level confidence interval 
for the 'true' required sample size}
\item{Nstar.low}{The lower endpoint of the desired level confidence interval 
for the 'true' required sample size}
\item{num.iters}{Number of iterations to compute. \code{num.iters = NULL} by 
default (i.e., all iterations are computed; users can specify 
\code{num.iters} = 1 for the first iteration.)}
\item{p}{The user-specified level of haplotype recovery. Default is 
\code{p} = 0.95.}
\item{perms}{The user-specified number of permutations (replications). Default 
is \code{perms} = 10000.}
\item{probs}{Haplotype frequency distribution vector}
\item{progress}{Should iteration results be outputted to the console? Default 
is TRUE.}
\item{prop.haps}{If \code{subset.haps} = TRUE, the user-specified proportion of haplotype labels to recover}
\item{prop.seqs}{If \code{subset.seqs} = TRUE, the user-specified proportion of 
DNA sequences to recover}
\item{ptm}{A timer to track progress of the algorithm in seconds}
\item{R}{The proportion of haplotypes recovered by the algorithm}
\item{R.low}{The lower endpoint of the desired level confidence interval for 
the 'true' fraction of haplotypes captured}
\item{R.up}{The upper endpoint of the desired level confidence interval for the 
'true' fraction of haplotypes captured}
\item{subset.haps}{Should a subsample of haplotype labels be taken? Default is 
FALSE.}
\item{subset.seqs}{Should a subsample of DNA sequences be taken? Default is 
FALSE.}
\item{X}{Mean number of specimens not sampled}
}

\examples{
# Returns the frequencies of each haplotype in the extrapolated sample 
max(envr$d$specs) * envr$probs

# Returns the extrapolated sample size corresponding to the dotted line 
# in the last iteration plot
envr$d[which(envr$d$means >= envr$p * envr$Hstar), ][1, 1]
}