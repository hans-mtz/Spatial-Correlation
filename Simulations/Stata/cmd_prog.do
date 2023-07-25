
quietly {
	program run_regs
		args x1
		mean `x1'
		scalar mean_x1 = r(table)[1,1]
		if mean_x1 == 1 {	
			reg y , robust
			gen p_val_HR = r(table)["pvalue","_cons"]
			gen t = r(table)["t","_cons"]
			scpc, cvs
			gen p_val_SCPC = e(scpcstats)["_cons","P>|t|"]
			gen t_SCPC = e(scpcstats)["_cons","t  "]
			gen t_SCPC_cv = e(scpccvs)["_cons","5%"]
		} 
		else {
			reg y `x1', robust
			gen p_val_HR = r(table)["pvalue","`x1'"]
			gen t = r(table)["t","`x1'"]
			scpc, cvs
			gen p_val_SCPC = e(scpcstats)["`x1'","P>|t|"]
			gen t_SCPC = e(scpcstats)["`x1'","t  "]
			gen t_SCPC_cv = e(scpccvs)["`x1'","5%"]
		}
	end

	run_regs x1
	
}

