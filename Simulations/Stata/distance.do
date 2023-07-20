clear all

import delimited "/Volumes/SSD Hans/Github/Spatial Correlation/Simulations/example.csv"


program distance
	args v1
	scalar c=75.167
// 	mkmat `1', matrix(vec)
	mata: getdistmat("`v1'")
	matrix Sigma = exp(-c*D)
// 	scalar n_rows = rowsof(vec)
// 	local n_r 250
// // 	display n_rows
// 	matrix d = J(n_rows,n_rows,0)
// // 	matrix list d
// 	forvalues i=1(1)`n_r' {
// 		forvalues j=1(1)`n_r' {
// // 			display `i'
// // 			matrix d[`i',`j'] = 1+0
// // 			matrix d[`i',`j'] = exp(-c*abs(vec[`i',...]-vec[`j',1]))
// 		}
// 	}
// 	matrix list d
	ereturn matrix Sigma
end
mata:

real matrix getdistmat(string scalar varname)
// computes matrix of distances from locations
// if external latlongflag=0, distances computed as norm, otherwise latitude / longitude and haversine formula for sphere with radius 1/Pi
{	
// 	external real scalar latlongflag
	
	real matrix mat
	real colvector s
	real scalar n,i //,c
// 	real matrix d
	st_view(s,.,varname)
	n=length(s)
	mat=J(n,n,0)
// 	if(latlongflag==0){
	for(i=1;i<=length(s);i++){
		mat[.,i]=sqrt(rowsum((s:-s[i,.]):^2))
	}
// 	}
// 	else
// 		{
// 		c=3.14159265359/180
// 		for(i=1;i<=n;i++){
// 			d=(.5*c)*(s:-s[i,.])
// 			mat[.,i]=asin(sqrt(sin(d[.,1]):^2  + cos(c*s[i,1])*(cos(c*s[.,1]):*(sin(d[.,2]):^2))))/3.14159265359
// 		}
// 	}	
// 	return(mat)
	st_matrix("D",mat)
}	
end

distance s_1
