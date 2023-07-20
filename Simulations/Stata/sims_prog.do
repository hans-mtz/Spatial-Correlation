clear all
program sims
	args n_s n_obs
	clear
// 	local n_s 5
	local r 0
	
	forvalues i=1(1)`n_s' {
		clear
		set obs `n_obs' // 250
		scalar beta = 1

		gen s_1 = runiform(0,1)
		gen e = rnormal(0,1)
		gen x1 = 1
		gen y = beta*x1 + e

		quietly reg y x, nocon robust
		quietly scpc
		local p_v=e(scpcstats)[1,4]
		if `p_v' > 0.05 {
			local ++r
		}
		local rej_freq = `r'/`n_s'
		
	}
	display "Rejection frequency was: `rej_freq'"
	
end

sims 100 250
