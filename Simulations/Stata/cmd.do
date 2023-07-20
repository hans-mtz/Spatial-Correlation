
quietly {
    set more off
    reg y x1, nocon robust
    gen p_val_HR = r(table)["pvalue","x1"]
	gen t = r(table)["t","x1"]
    scpc, cvs
    gen p_val_SCPC = e(scpcstats)["x1","P>|t|"]
	gen t_SCPC = e(scpcstats)["x1","t  "]
	gen t_SCPC_cv = e(scpccvs)["x1","5%"]
}



