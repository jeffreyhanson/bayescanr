#' @include bayescanr-internal.R misc.R generics.R
NULL

#' BayeScanData: An S4 class to represent input data for BayeScan
#'
#' This class stores input data for the BayeScan program.
#'
#' @slot matrix \code{matrix} with binary loci data. Each row is for a different sample; each column is for a different loci.
#' @slot primers \code{character} names of the selective primers used for the AFLPs followed by the location on the gel.
#' @slot populations \code{character} name of populations that each sample is from.
#' @slot labels \code{character} labels for each sample.
#' @seealso \code{\link{BayeScanData}}.
#' @export
setClass(
	"BayeScanData",
	representation(
		matrix='matrix',
		primers='character',
		populations='character',
		labels='character'
	),
	validity=function(object) {
		# matrix
		expect_true(all(object@matrix[] %in% c(0, 1, NA)))
		# check that dimensions > 0
		expect_true(all(dim(object@matrix)>0))
		# primers
		expect_is(object@primers,'character')
		expect_true(all(!is.na(object@primers)))
		# populations
		expect_is(object@populations,'character')
		expect_true(all(!is.na(object@populations)))
		# labels
		expect_is(object@labels,'character')
		expect_true(all(!is.na(object@labels)))
		# cross-object checks
		expect_equal(length(object@primers),ncol(object@matrix))
		expect_equal(length(object@populations),nrow(object@matrix))
		expect_equal(length(object@labels),nrow(object@matrix))
		return(TRUE)
	}
)

#' Create BayeScanData object
#'
#' This function creates a new \code{BayeScanData} object.
#'
#' @param matrix \code{matrix} with binary loci data. Each row is for a different sample; each column is for a different loci.
#' @param primers \code{character} names of the selective primers used for the AFLPs followed by the location on the gel.
#' @param populations \code{character} name of populations that each sample is from.
#' @param labels \code{character} labels for each sample.
#' @seealso \code{\link{BayeScanData-class}}.
#' @export
BayeScanData<-function(matrix, primers, populations, labels) {
	x<-new("BayeScanData", matrix=matrix, primers=primers, populations=populations, labels=labels)
	validObject(x, test=FALSE)
	return(x)
}

#' @rdname n.loci
#' @method n.loci BayeScanData
#' @export
n.loci.BayeScanData <- function(x) {
	return(ncol(x@matrix))
}

#' @rdname n.pop
#' @method n.pop BayeScanData
#' @export
n.pop.BayeScanData <- function(x) {
	return(length(unique(x@populations)))
}

#' @rdname n.samples
#' @method n.samples BayeScanData
#' @export
n.samples.BayeScanData <- function(x) {
	return(nrow(x@matrix))
}

#' @rdname pop.names
#' @method pop.names BayeScanData
#' @export
pop.names.BayeScanData <- function(x) {
	return(unique(x@populations))
}

#' @rdname sample.pops
#' @method sample.pops BayeScanData
#' @export
sample.pops.BayeScanData <- function(x) {
	return(x@populations)
}

#' @rdname sample.pops
#' @method sample.pops<- BayeScanData
#' @export
`sample.pops<-.BayeScanData` <- function(x, value) {
	x@populations <- value
	return(x)
}

#' @rdname sample.labels
#' @method sample.labels BayeScanData
#' @export
sample.labels.BayeScanData <- function(x) {
	return(x@labels)
}

#' @rdname sample.labels
#' @method sample.labels<- BayeScanData
#' @export
`sample.labels<-.BayeScanData` <- function(x, value) {
	x@labels <- value
	return(x)
}

#' @rdname pop.subset
#' @method pop.subset BayeScanData
#' @export
pop.subset.BayeScanData <- function(x, populations) {
	pos <- which(x@populations %in% populations)
	return(
		BayeScanData(
			matrix=x@matrix[pos,,drop=FALSE],
			populations=x@populations[pos],
			primers=x@primers,
			labels=x@labels[pos]
		)
	)
}

#' @rdname loci.subset
#' @method loci.subset BayeScanData
#' @export
loci.subset.BayeScanData <- function(x, loci) {
	if (is.character(loci)) {
		stop('argument to loci must be numeric when x is BayeScanData')
	} else {
		return(
			BayeScanData(
				matrix=x@matrix[,loci,drop=FALSE],
				populations=x@populations,
				primers=x@primers[loci],
				labels=x@labels
			)
		)
	}
}

#' @rdname sample.subset
#' @method sample.subset BayeScanData
#' @export
sample.subset.BayeScanData <- function(x, samples) {
	if (is.character(samples))
		samples <- match(samples, x@labels)
	return(
		BayeScanData(
			matrix=x@matrix[samples,,drop=FALSE],
			populations=x@populations[samples],
			primers=x@primers,
			labels=x@labels[samples]
		)
	)
}


#' Read FSTAT data for BayeScan
#'
#' This function reads FSTAT data for BayeScan.
#'
#' @param x \code{character} file path with data
#' @seealso \code{BayeScanData}.
#' @return \code{BayeScanData}.
#' @examples
#' x <- read.BayeScanData(system.file('extdata', 'example_fstat_aflp.dat', package='bayescanr'))
#' @export
read.BayeScanData <- function(x) {
	# read first line
	meta <- as.numeric(strsplit(readLines(x, n=1), '\t')[[1]])
	# read in primer labels
	primers <- scan(x, skip=1, n=meta[2], what='character', quiet=TRUE)
	# read in population ids
	dat <- fread(x, skip=meta[2]+1, data.table=FALSE)
	na.col <- apply(dat, 2, function(x) {all(is.na(x))})
	if (any(na.col))
		dat <- dat[,-which(na.col),drop=FALSE]
	# substitute 0, 1, NA
	mat <- as.matrix(dat[,-1,drop=FALSE])
	mat[which(mat[]==0)] <- NA
	mat[which(mat[]==1)] <- 0
	mat[which(mat[]==2)] <- 1
	# return object
	return(
		BayeScanData(matrix=mat, primers=primers, populations=as.character(dat[[1]]), labels=as.character(seq_len(nrow(mat))))
	)
}

#' Write data for BayeScan
#'
#' This function writes data for BayeScan.
#'
#' @param x \code{BayeScanData} object.
#' @param file \code{character} file path to write data.
#' @seealso \code{BayeScanData}.
#' @details This code is based on the R code distributed with BayeScan (version 2.1) for converting binary FSTAT files to the BayeScan data format.
#' @examples
#' x <- read.BayeScanData(system.file('extdata', 'example_fstat_aflp.dat', package='bayescanr'))
#' write.BayeScanData(x, tempfile(fileext='txt'))
#' @export
write.BayeScanData <- function(x,file) {
	# init
	out <- c(
		'[loci]=',n.loci(x),'\n\n',
		'[populations]=',n.pop(x),'\n\n'
	)
	# main
	loci.ids <- paste0(seq_len(n.loci(x)), ' ')
	curr.pops.names <- pop.names(x)
	for (i in seq_len(n.pop(x))) {
		# init
		out <- c(out,"[pop]=",i,'\n')
		curr.subset <- pop.subset(x, curr.pops.names[i])
		n.with <- paste0(colSums(curr.subset@matrix, na.rm=TRUE), '\n')
		n.total <- paste0(colSums(!is.na(curr.subset@matrix)), ' ')
		# append data
		out <- c(
			out,
			c(t(matrix(c(loci.ids,n.total,n.with),ncol=3))),
			'\n'
		)
	}
	# exports
	cat(out, file=file, sep='')
	return(invisible())
}

#' @method nmds BayeScanData
#' @rdname nmds
#' @export
nmds.BayeScanData <- function(x, max.stress=0.1, min.k=2, max.k=Inf, metric='gower', type='all', ...) {
	# init
	stopifnot(min.k<=max.k)
	if (type!='all')
		stop('argument to type must be all when x inherits from BayeScanData')
	# prelim
	dist.mtx <- daisy(
		cbind(as.data.frame(x@matrix==1),1),
		metric=metric,
		type=list(asymm=seq_len(n.loci(x)))
	)
	# main
	curr.stress <- Inf
	curr.k <- min.k
	# find nmds with suitable k
	while (curr.stress > max.stress & curr.k <= max.k) {
		print('here')
		curr.nmds <- metaMDS(
			comm=dist.mtx,
			k=curr.k,
			wascores=FALSE,
			autotransform=FALSE,
			noshare=FALSE,
			...
		)
		curr.stress <- curr.nmds$stress
		curr.k <- curr.k + 1
	}
	# reutrn object
	return(curr.nmds)
}

#' @method print BayeScanData
#' @rdname print
#' @export
print.BayeScanData <- function(x, ..., header=TRUE) {
	if (header)
		cat("BayeScanData object.\n")
	cat('  populations:',nrow(x@matrix),'\n')
	cat('  loci:',ncol(x@matrix),'\n')
}

#' @rdname show
#' @method show BayeScanData
#' @export
setMethod(
	'show',
	'BayeScanData',
	function(object)
		print.BayeScanData(object)
)
