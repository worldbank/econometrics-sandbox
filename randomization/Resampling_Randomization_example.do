**** Example code for the resampling techniques to test for randomization used  in Does Professor Quality Matter? Evidence from Random Assignment of Students to Professors by Scott E. Carrell, and James E. West
  
****Initialize a matrix to store p-values
matrix A=0,0,0
**** set a seed to get the same results each time you run the code.
set seed 123456789
***** In the following example, the randomization was stratified by year, so we do the resampling and checking for each year.
forvalues year = 2004(1)2016 {
   display `year'
   preserve
**** Keeps only the relevant year for the loop   
   quietly drop if year != `year'
***** drops observations for which the covariate of interest is missing.    
   quietly drop if covariate ==.
   gen select = 0
   gen udraw = 0
	quietly levelsof group_id, local(levels)
	foreach x of local levels {
	   display `x'
 	   quietly replace select = 0
	   quietly replace select = 1 if group_id== `x'
	   ***** obtain the mean or sum of the covariate of interest for the observed group
	   quietly summarize covariate if select
	   scalar nsec = r(N)
	   scalar sumorig = r(sum)
	   scalar accum = 0
	   ***** Redraw group 10000 times, each time recording whether the mean or sum of covariate of interest was lower than for the observed group
	   forvalues i = 1(1)10000 {
		quietly replace udraw = uniform()
		sort udraw, stable
		quietly summarize covariate if _n <= nsec
		scalar sumfake=r(sum)
		scalar accum = accum + (sumfake < sumorig)  
		}
		***** Save the year, group, p-value (number of times lower/total resamples)
		matrix A=A\ `year' , `x',(accum/10000)
	}
	restore
}

svmat A
keep A1 A2 A3 
drop if _n==1
rename A1 year
rename A2 group_id
rename A3 pvalue

****The program below runs a chi-square and ksmirnow test for each year
capture program drop KsChiSq
program KsChiSq
 forvalues year = 2004(1)2016 {
     preserve
	quietly drop if year != `year'
	quietly drop if pvalue == 0 | pvalue == 1 
		  scalar teststat = 0
		  quietly sum pvalue
		  scalar nobs = r(N)
		  forvalues i = 0.1(0.1)1 {
		    quietly sum pvalue if ((pvalue >= (`i'-0.1)) & (pvalue < `i'))
		    scalar teststat = teststat + ((r(N)-0.1*nobs)^2)/(nobs*0.1)
		   }
		local pval = chi2(10,teststat)
		capture quietly ksmirnov pvalue = pvalue
		local pt = r(p_cor)
		  display "`year', `pt', `pval'"
		restore
 }
end
}
************* For each year compare p-values to 0.05 to determine if the tests failed. The first column is year, the second is the p-value from the chi-squared test, the third is the p-value from the ksmirnov test.
KsChiSq
