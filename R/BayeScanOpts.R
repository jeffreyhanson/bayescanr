#' @include bayescanr-internal.R misc.R generics.R
NULL

#' BayeScanOpts: An S4 class to represent BayeScan parameters
#'
#' This class stores input parameters for the BayeScan program.
#'
#' @slot threads \code{numeric} number of threads used for computation. Defaults to 1.
#' @slot reps \code{numeric} number of replicate runs. Defaults to 3.
#' @slot n \code{numeric} number of outputted iterations. Defaults to 5000.
#' @slot thin \code{numeric} thinning interval size. Defaults to 10.
#' @slot nbp \code{numeric} number of pilot runs. Defaults to 20.
#' @slot pilot \code{numeric} length of pilot runs. Defaults to 5000.
#' @slot burn \code{numeric} burn-in length. Defaults to 50000.
#' @slot fdr \code{numeric} false discovery rate used to classify loic as adaptive. Defaults to 0.1.
#' @seealso \code{\link{BayeScanOpts}}.
#' @export
setClass(
	"BayeScanOpts",
	representation(
		threads="numeric",
		reps="numeric",
		n="numeric",
		thin="numeric",
		nbp="numeric",
		pilot="numeric",
		burn="numeric",
		fdr="numeric"
	),
	prototype=list(
		threads=1,
		reps=3,
		n=5000,
		thin=10,
		nbp=20,
		pilot=5000,
		burn=50000,
		fdr=0.1
	),
	validity=function(object) {
		# check that parameters are greater than zero and not NA
		sapply(
			slotNames(object),
			function(x) {
				expect_true(is.finite(slot(object, x)))
				expect_true(slot(object, x)>0)
				return(invisible())
		})
		# check that threads is not higher than the number of threads on cpu
		expect_true(object@threads <= detectCores())
		return(TRUE)
	}
)

#' Create BayeScanOpts object
#'
#' This function creates a new \code{BayeScanOpts} object.
#'
#' @param threads \code{numeric} number of threads used for computation. Defaults to 1.
#' @param reps \code{numeric} number of replicate runs. Defaults to 3.
#' @param n \code{numeric} number of outputted iterations. Defaults to 5000.
#' @param thin \code{numeric} thinning interval size. Defaults to 10.
#' @param nbp \code{numeric} number of pilot runs. Defaults to 20.
#' @param pilot \code{numeric} length of pilot runs. Defaults to 5000.
#' @param burn \code{numeric} burn-in length. Defaults to 50000.
#' @param fdr \code{numeric} false discovery rate used to classify loic as adaptive. Defaults to 0.1.
#' @seealso \code{\link{BayeScanOpts-class}}.
#' @export
BayeScanOpts<-function(threads=1, reps=3, n=5000, thin=10, nbp=20, pilot=5000, burn=50000, fdr=0.1) {
	x<-new("BayeScanOpts", threads=threads, reps=reps, n=n, thin=thin, nbp=nbp, pilot=pilot, burn=burn, fdr=fdr)
	validObject(x, test=FALSE)
	return(x)
}

#' @method print BayeScanOpts
#' @rdname print
#' @export
print.BayeScanOpts=function(x, ..., header=TRUE) {
	if (header)
		cat("BayeScanOpts object.\n")
	cat('  threads:',x@threads,'\n')
	cat('  reps:',x@reps,'\n')
	cat('  n:',x@n,'\n')
	cat('  thin:',x@thin,'\n')
	cat('  nbp:',x@nbp,'\n')
	cat('  pilot:',x@pilot,'\n')
	cat('  burn:',x@burn,'\n')
	cat('  fdr:',x@fdr,'\n')
}

#' @rdname show
#' @export
setMethod(
	'show',
	'BayeScanOpts',
	function(object)
		print.BayeScanOpts(object)
)
