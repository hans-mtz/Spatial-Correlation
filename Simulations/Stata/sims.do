// Parameters
clear all
set obs 250

scalar n_s = 100
scalar beta = 1

// DGP
gen s_1 = runiform(0,1)
gen e = rnormal(0,1)
gen x1 = 1
gen y = beta*x1 + e

// Regression
reg y x1, nocon robust

// Spatial correlation comand
scpc

// retrieving the p-value
scalar p_v=e(scpcstats)[1,4]

// To see all values stored in e() 
ereturn list all

// Simulation program 



program testing
	clear
	forvalues i=1(1)`1' {
		display " Printing `i' of `1'"
		
	}
	
end

testing 9
