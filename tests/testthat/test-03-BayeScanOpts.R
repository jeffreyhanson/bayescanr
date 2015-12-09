test_that('BayeScanOpts', {
	# tests implicit
	x <- BayeScanOpts(threads=1, reps=3, n=5000, thin=10, nbp=20, pilot=5000, burn=50000, fdr=0.1)
	# methods
	print(x)
	x
})

