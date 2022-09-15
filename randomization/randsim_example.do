/*******************************************************************************
This do file is created for a blog "We randomized! But did we really, though?"
It provides three example checks to test whether a randomization process is 
working as intended.

randsim_example.do

2022.09.14 Created
*******************************************************************************/
set seed 155159 // chosen from random.org Min: 1, Max: 1000000 2022-08-11 03:47:57 UTC 

forvalues i=1/200{
	sysuse auto, clear
	gen sim = `i'
	randtreat, gen(T) strata(foreign) misfits(global)
	reg mpg T i.foreign, vce(hc3)
	gen coef = _b[T] in 1
	test T
	gen p=r(p) in 1
	cap append using `dat'
	tempfile dat
	save `dat'
}	

cap graph drop _all

// (1) Check if the distribution of the saved p-values is uniform
kdensity p, name(p) graphregion(color(white))

// (2) Check if the distribution of the saved coefficients is mean zero
kdensity coef, name(coef) graphregion(color(white))

// (3) Check if any unit has a 0 probability of being assigned to T across draws
collapse T, by(make)
tab make if T==0
