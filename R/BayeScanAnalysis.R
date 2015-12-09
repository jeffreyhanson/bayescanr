#' @include bayescanr-internal.R misc.R generics.R BayeScanOpts.R BayeScanData.R BayeScanResults.R
NULL

#' BayeScanAnalysis: An S4 class to represent inputs and outputs from BayeScan
#'
#' This class stores input data and associated output results from the BayeScan program.
#'
#' @slot opts \code{BayeScanOpts} object with parameters used to run BayeScan .
#' @slot data \code{BayeScanData} object with input data used for analysis.
#' @slot results \code{BayeScanResults} object with results from analysis.
#' @seealso \code{\link{BayeScanAnalysis}}.
#' @export
setClass(
	"BayeScanAnalysis",
	representation(
		opts='BayeScanOpts',
		data='BayeScanData',
		results='BayeScanResults'
	),
	validity=function(object) {
		## opts
		# checks are internal
		## data
		# checks are internal
		## results
		# checks are internal
		## cross-object checks
		expect_equal(n.pop(object@data), n.pop(object@results))
		expect_equal(n.loci(object@data), n.loci(object@results))
		return(TRUE)
	}
)

#' Create BayeScanAnalysis object
#'
#' This function creates a new \code{BayeScanAnalysis} object.
#'
#' @param opts \code{BayeScanOpts} object with parameters used to run BayeScan .
#' @param data \code{BayeScanData} object with input data used for analysis.
#' @param results \code{BayeScanResults} object with results from analysis.
#' @seealso \code{\link{BayeScanAnalysis-class}}, \code{\link{BayeScanData}}, \code{\link{BayeScanData}}, \code{\link{BayeScanResults}}.
#' @export
BayeScanAnalysis<-function(opts, data, results) {
	x<-new("BayeScanAnalysis", opts=opts, data=data, results=results)
	validObject(x, test=FALSE)
	return(x)
}

#' @rdname n.loci
#' @method n.loci BayeScanAnalysis
#' @export
n.loci.BayeScanAnalysis <- function(x) {
	return(n.loci(x@data))
}

#' @rdname n.pop
#' @method n.pop BayeScanAnalysis
#' @export
n.pop.BayeScanAnalysis <- function(x) {
	return(n.pop(x@data))
}

#' @rdname n.samples
#' @method n.samples BayeScanAnalysis
#' @export
n.samples.BayeScanAnalysis <- function(x) {
	return(n.samples(x@data))
}

#' @rdname pop.names
#' @method pop.names BayeScanAnalysis
#' @export
pop.names.BayeScanAnalysis <- function(x) {
	return(pop.names(x@data))
}

#' @rdname sample.pops
#' @method sample.pops BayeScanAnalysis
#' @export
sample.pops.BayeScanAnalysis <- function(x) {
	return(sample.pops(x@data))
}

#' @rdname loci.subset
#' @method loci.subset BayeScanAnalysis
#' @export
loci.subset.BayeScanAnalysis <- function(x, loci) {
	if (is.character(loci)) {
		if (loci=='all') {
			return(x@data)
		} else {
			return(loci.subset(x@data, x@results@summary$type==loci))
		}
	} else {
		return(loci.subset(x@data, loci))
	}
}


#' Run BayeScan
#'
#' This function analysis data using BayeScan.
#'
#' @param x \code{BayeScanData} object.
#' @inheritParams BayeScanOpts
#' @param dir \code{character} with directory to use for analysis.
#' @param clean \code{logical} should input and output files be deleted after analysis is finished?
#' @seealso \code{BayeScanData}, \code{BayeScanOpts}.
#' @examples
#' # run BayeScan using low number of iterations
#' dat <- read.BayeScanData(system.file('extdata', 'example_fstat_aflp.dat', package='bayescanr'))
#' x <- run.BayeScan(dat, threads=1, n=50, thin=1, nbp=10, pilot=10, burn=10)
#' @export
run.BayeScan<-function(x, threads=1, reps=3, n=5000, thin=10, nbp=20, pilot=5000, burn=50000, fdr=0.1, dir=tempdir(), clean=TRUE) {
	## initialization
	# argument checks
	opts <- BayeScanOpts(threads=threads, reps=reps, n=n, thin=thin, nbp=nbp, pilot=pilot, burn=burn, fdr=fdr)
	expect_is(x, 'BayeScanData')
	# set BayeScan file path
	bayescan.path <- switch(
		Sys.info()['sysname'],
		'Linux'=system.file('bin', 'BayeScan2.1_linux64bits', package='bayescanr'),
		'Darwin'=system.file('bin', 'bayescan_2.1_macos64bits', package='bayescanr'),
		'Windows'=system.file('bin', 'bayescan_2.1_win32bits_cmd_line.exe', package='bayescanr'),
	)
	# update permissions
	if (!grepl(basename(bayescan.path), 'win'))
		system(paste0('chmod 700 ',bayescan.path))
	### main processing
	# write data to file
	dat.path <- tempfile(tmpdir=dir, fileext='.txt')
	write.BayeScanData(x, dat.path)
	# run BayesScan analysis
	replicates <- lapply(seq_len(opts@reps), function(i) {
		system(
			paste0(
				bayescan.path, ' ',
				dat.path,
				' -od ',dir,
				' -threads ',opts@threads,
				' -n ',opts@n,
				' -thin ',opts@thin,
				' -nbp ',opts@nbp,
				' -pilot ',opts@pilot,
				' -burn ',opts@burn
			)
		)
		return(read.BayeScanReplicate(dat.path,dir,fdr=opts@fdr))
	})
	## exports
	# construct BayeScanAnalysis object
	return(
		BayeScanAnalysis(
			opts=opts,
			data=x,
			results=BayeScanResults(replicates=replicates)
		)
	)
}

#' @method mds BayeScanAnalysis
#' @rdname mds
#' @export
mds.BayeScanAnalysis <- function(x, metric='gower', type='all', ...) {
	return(
		mds.BayeScanData(
			loci.subset(x, type),
			metric,
			...
		)
	)
}

#' @method print BayeScanAnalysis
#' @rdname print
#' @export
print.BayeScanAnalysis=function(x, ..., header=TRUE) {
	if (header)
		cat("BayeScanAnalysis object.\n\n")
	cat('Options','\n')
	print(x@data, header=FALSE)
	cat('Data','\n')
	print(x@data, header=FALSE)
	cat('Results','\n')
	print(x@results, header=FALSE)
}

#' @rdname show
#' @export
setMethod(
	'show',
	'BayeScanAnalysis',
	function(object)
		print.BayeScanAnalysis(object)
)

