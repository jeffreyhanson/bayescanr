context('misc functions')

test_that('compatibility with coda package', {
	library(coda)
	x <- mcmc.list(mcmc(rnorm(100)),mcmc(rnorm(100)))
	g <- gelman.diag(x)
	print(g)
})
