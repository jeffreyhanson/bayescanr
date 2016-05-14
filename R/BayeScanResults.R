#' @include bayescanr-internal.R misc.R generics.R BayeScanReplicate.R
NULL

#' BayeScanResults: An S4 class to results from BayeScan
#'
#' This class stores results from the BayeScan program.
#'
#' @slot summary \code{data.frame} object containing overall results from BayeScan replicates.
#' @slot replicates \code{list} of \code{BayeScanReplicate} objects.
#' @seealso \code{\link{BayeScanResults}}, \code{\link{BayeScanReplicate}}.
#' @export
setClass(
	"BayeScanResults",
	representation(
		summary='data.frame',
		replicates='list'
	),
	validity=function(object) {
		# check that all replicates have the same properties
		expect_equal(length(unique(sapply(object@replicates, n.loci))), 1)
		expect_equal(length(unique(sapply(object@replicates, n.pop))), 1)
		return(TRUE)
	}
)

#' Create BayeScanResults object
#'
#' This function creates a new \code{BayeScanResults} object.
#'
#' @param summary \code{data.frame} object containing overall results from BayeScan replicates.
#' @param replicates \code{list} of \code{BayeScanReplicate} objects.
#' @seealso \code{\link{BayeScanReplicate-class}}.
#' @return \code{\link{BayeScanReplicate}}.
#' @export
BayeScanResults<-function(summary=NULL, replicates) {
	# compute summary results
	if (is.null(summary)) {
		if (length(replicates)>1) {
			summary <- data.frame(
				loci=replicates[[1]]@fst$loci,
				mean_prob=rowMeans(sapply(replicates, function(x) {x@fst$prob})),
				mean_log10_PO=rowMeans(sapply(replicates, function(x) {x@fst$log10_PO})),
				mean_qval=rowMeans(sapply(replicates, function(x) {x@fst$qval})),
				mean_alpha=rowMeans(sapply(replicates, function(x) {x@fst$alpha})),
				mean_fst=rowMeans(sapply(replicates, function(x) {x@fst$fst})),
				type=apply(
					as.matrix(sapply(replicates, function(x) {as.character(x@fst$type)})),
					1,
					function(x) {c('neutral','adaptive')[1+all(x=='adaptive')]}
				)
			)
		} else {
			summary <- replicates[[1]]@fst
		}
	}
	# return new object
	x<-new("BayeScanResults", summary=summary, replicates=replicates)
	validObject(x, test=FALSE)
	return(x)
}

#' @rdname n.loci
#' @method n.loci BayeScanResults
#' @export
n.loci.BayeScanResults <- function(x) {
	return(nrow(x@summary))
}

#' @rdname n.pop
#' @method n.pop BayeScanResults
#' @export
n.pop.BayeScanResults <- function(x) {
	return((ncol(x@replicates[[1]]@mcmc)-2)/2)
}

#' @rdname n.samples
#' @method n.samples BayeScanResults
#' @export
n.samples.BayeScanResults <- function(x) {
	stop('BayeScanResults does not store this information.')
	return(invisible())
}

#' @rdname pop.names
#' @method pop.names BayeScanResults
#' @export
pop.names.BayeScanResults <- function(x) {
	stop('BayeScanResults does not store population names.')
	return(invisible())
}

#' @rdname sample.pops
#' @method sample.pops BayeScanResults
#' @export
sample.pops.BayeScanResults <- function(x) {
	stop('BayeScanResults does not store population names.')
	return(invisible())
}

#' @method print BayeScanResults
#' @rdname print
#' @export
print.BayeScanResults=function(x, ..., header=TRUE) {
	if (header)
		cat("BayeScanResults object.\n")
	cat('  adaptive loci:',sum(x@summary$type=='adaptive'),'\n')
	cat('  neutral loci:',sum(x@summary$type=='neutral'),'\n')
	cat('  Gelman-Ruben R:', gelman.diag(x)$psrf[[1]], '\n')
}

#' @method traceplot BayeScanResults
#' @rdname traceplot
#' @export
traceplot.BayeScanResults <- function(x, ...) {
	# extract logliks
	ll <- data.frame(iteration=x@replicates[[1]]@mcmc$iteration)
	ll <- cbind(ll, data.frame(sapply(x@replicates, function(y) {y@mcmc$logL})))
	ll <- gather(ll, chain, loglik, -iteration)
	ll$chain <- gsub('X', '', ll$chain, fixed=TRUE)
	# make plot
	ggplot(data=ll, aes(x=iteration, y=loglik, color=chain)) +
		geom_line() + xlab('Iteration') + ylab('Negative loglikelihood')
}

#' @method gelman.diag BayeScanResults
#' @rdname gelman.diag
#' @export
gelman.diag.BayeScanResults <- function(x, ...) {
	if(length(x@replicates)>1)
		return(
			coda::gelman.diag(
				do.call(
					mcmc.list,
					lapply(x@replicates, function(y) {mcmc(y@mcmc$logL, thin=diff(y@mcmc$iteration[2:1]))})
				)
			)
		)
	# return gelman.diag object with NAs if only one chain
	return(structure(list(psrf = structure(c(NA_real_, NA_real_), .Dim = 1:2, .Dimnames = list(NULL, c("Point est.", "Upper C.I."))),
		mpsrf = NULL), .Names = c("psrf", "mpsrf"), class = "gelman.diag"))
}

#' @rdname show
#' @export
setMethod(
	'show',
	'BayeScanResults',
	function(object)
		print.BayeScanResults(object)
)

 
