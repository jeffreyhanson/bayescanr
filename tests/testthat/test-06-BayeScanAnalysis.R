context('BayeScanAnalysis')

test_that('run.BayeScan', {
	# make BayeScan object
	bd <- read.BayeScanData(system.file('extdata', 'example_fstat_aflp.dat', package='bayescanr'))
	bs <- run.BayeScan(bd, fdr=0.5,threads=1,n=10,thin=1,nbp=5,pilot=5,burn=50,reps=2)
	# methods
	nmds(bs, metric='gower', type='neutral', trymax=2, max.stress=0.3)
	print(bs)
	bs
	n.loci(bs)
	n.pop(bs)
	gelman.diag(bs)
	traceplot(bs)
})

