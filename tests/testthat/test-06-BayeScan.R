test_that('run.BayeScan', {
	# make BayeScan object
	bd <- read.BayeScanData(system.file('extdata', 'example_fstat_aflp.dat', package='bayescanr'))
	bs <- run.BayeScan(bd, threshold=0.5,threads=1,n=10,thin=1,nbp=5,pilot=5,burn=50)
	# methods
	mds(bs, metric='gower', k=2, type='neutral', trymax=2)
	print(bs)
	bs
})

