test_that('read.BayeScanResults', {
	# create data
	dir<-tempdir()
	path<-tempfile(tmpdir=dir, fileext='.txt')
	bo <-BayeScanOpts(fdr=0.5,threads=1,reps=1,n=10,thin=1,nbp=5,pilot=5,burn=50)
	bd <- read.BayeScanData(system.file('extdata', 'example_fstat_aflp.dat', package='bayescanr'))
	write.BayeScanData(bd,path)
	# identify bayescan path
	bayescan.path <- switch(
		Sys.info()['sysname'],
		'Linux'=system.file('bin', 'BayeScan2.1_linux64bits', package='bayescanr'),
		'Darwin'=system.file('bin', 'BayeScan2.1_macos64bits', package='bayescanr'),
		'Windows'=system.file('bin', 'BayeScan2.1_win32bits_cmd_line.exe', package='bayescanr')
	)
	# update permissions
	if (!grepl(basename(bayescan.path), 'win'))
		system(paste0('chmod 700 ',bayescan.path))
	# run BayeScan
	system(
		paste0(
			bayescan.path, ' ',
			path,
			' -od ',dir,
			' -threads ',bo@threads,
			' -n ',bo@n,
			' -thin ',bo@thin,
			' -nbp ',bo@nbp,
			' -pilot ',bo@pilot,
			' -burn ',bo@burn
		)
	)
	# try reading results back into R
	results <- BayeScanResults(replicates=list(read.BayeScanReplicate(path,dir,bo@fdr)))
	# methods
	print(results)
	results
	n.loci(results)
	n.pop(results)
})
