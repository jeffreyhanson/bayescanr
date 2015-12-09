#' @include bayescanr-internal.R misc.R generics.R
NULL

#' BayeScanReplicate: An S4 class to results from BayeScan
#'
#' This class stores results from the BayeScan program.
#'
#' @slot fst \code{data.frame} object containing output results ('filename_fst.txt').
#' @slot mcmc \code{data.frame} object containing Markov chain Monte Carlo information ('filename.sel').
#' @slot acceptance.rate \code{data.frame} object with information on the evolution of the acceptance rate ('filename_AccRte.txt').
#' @slot verification \code{character} object with verification diagnostics ('file_Verif.txt').
#' @seealso \code{\link{BayeScanReplicate}}.
#' @export
setClass(
	"BayeScanReplicate",
	representation(
		fst='data.frame',
		mcmc='data.frame',
		acceptance.rate='data.frame',
		verification='character'
	),
	validity=function(object) {
		# fst
		expect_is(object@fst, 'data.frame')
		expect_equal(names(object@fst), c("loci", "prob", "log10_PO", "qval", "alpha", "fst", "type"))
		for (i in c("loci", "prob", "log10_PO", "qval", "alpha", "fst"))
			expect_is(object@fst[[i]], c('numeric','integer'))
		expect_is(object@fst[['type']], 'factor')
		# mcmc
		expect_is(object@mcmc, 'data.frame')
		for (i in names(object@mcmc))
			expect_is(object@mcmc[[i]], c('numeric','integer'))
		# acceptance.rate
		expect_is(object@acceptance.rate, 'data.frame')
		expect_equal(names(object@acceptance.rate), c("beta", "ances", "freq"))
		for (i in names(object@acceptance.rate))
			expect_is(object@acceptance.rate[[i]], c('numeric','integer'))
		# verification
		expect_is(object@verification, 'character')
		expect_equal(length(object@verification), 1L)
		return(TRUE)
	}
)

#' Create BayeScanReplicate object
#'
#' This function creates a new \code{BayeScanReplicate} object.
#'
#' @param fst \code{data.frame} object containing output results ('filename_fst.txt').
#' @param mcmc \code{data.frame} object containing Markov chain Monte Carlo information ('filename.sel').
#' @param acceptance.rate \code{data.frame} object with information on the evolution of the acceptance rate ('filename_AccRte.txt').
#' @param verification \code{character} object with verification diagnostics ('file_Verif.txt').
#' @seealso \code{\link{BayeScanReplicate-class}}.
#' @return \code{\link{BayeScanReplicate}}.
#' @export
BayeScanReplicate<-function(fst, mcmc, acceptance.rate, verification) {
	x<-new("BayeScanReplicate", fst=fst, mcmc=mcmc, acceptance.rate=acceptance.rate, verification=verification)
	validObject(x, test=FALSE)
	return(x)
}

#' @rdname n.loci
#' @method n.loci BayeScanReplicate
#' @export
n.loci.BayeScanReplicate <- function(x) {
	return(nrow(x@fst))
}

#' @rdname n.pop
#' @method n.pop BayeScanReplicate
#' @export
n.pop.BayeScanReplicate <- function(x) {
	return((ncol(x@mcmc)-2)/2)
}

#' @rdname n.samples
#' @method n.samples BayeScanReplicate
#' @export
n.samples.BayeScanReplicate <- function(x) {
	stop('BayeScanReplicate does not store this information.')
	return(invisible())
}

#' @rdname pop.names
#' @method pop.names BayeScanReplicate
#' @export
pop.names.BayeScanReplicate <- function(x) {
	stop('BayeScanReplicate does not store population names.')
	return(invisible())
}

#' @rdname sample.pops
#' @method sample.pops BayeScanReplicate
#' @export
sample.pops.BayeScanReplicate <- function(x) {
	stop('BayeScanReplicate does not store population names.')
	return(invisible())
}

#' Read BayeScan run
#'
#' This function reads the results of a single run of the BayeScan program.
#'
#' @param file \code{character} file path of input file.
#' @param dir \code{dir} directory with output files.
#' @param fdr \code{numeric} false discovery rate threshold to classify loci as adaptive.
#' @seealso \code{\link{BayeScanReplicate-class}}.
#' @return \code{\link{BayeScanReplicate}}.
#' @export
read.BayeScanReplicate<-function(file, dir, fdr=0.95) {
	# avoid cran note
	prob <- NULL
	# return object
	return(
		BayeScanReplicate(
			fst=base::transform(
				`names<-`(
					fread(gsub('\\.txt', '_fst.txt', file.path(dir, basename(file))), skip=1, data.table=FALSE),
					c('loci','prob','log10_PO','qval','alpha','fst')
				),
				type=c('neutral','adaptive')[(qval<fdr)+1]
			),
			mcmc=`names<-`(
				fread(gsub('\\.txt', '.sel', file.path(dir, basename(file))), skip=1, data.table=FALSE),
				c('iteration',strsplit(readLines(gsub('\\.txt', '.sel', file.path(dir, basename(file))),n=1),' ')[[1]])
			),
			acceptance.rate=fread(gsub('\\.txt', '_AccRte.txt', file.path(dir, basename(file))), data.table=FALSE),
			verification=paste(readLines(gsub('\\.txt', '_Verif.txt', file.path(dir, basename(file)))), collapse='\n')
		)
	)
}

#' @method print BayeScanReplicate
#' @rdname print
#' @export
print.BayeScanReplicate=function(x, ..., header=TRUE) {
	if (header)
		cat("BayeScanReplicate object.\n")
	cat('  adaptive loci:',sum(x@fst$type=='adaptive'),'\n')
	cat('  neutral loci:',sum(x@fst$type=='neutral'),'\n')
}

#' @rdname show
#' @export
setMethod(
	'show',
	'BayeScanReplicate',
	function(object)
		print.BayeScanReplicate(object)
)

 
