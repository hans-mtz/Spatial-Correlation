clear all
sysuse auto, clear
gen s_1=rnormal(0,1)
gen s_2=rnormal(0,1)
reg mpg weight length, robust
scpc 

matrix list e(scpccvs)
matrix list e(scpcstats)

matrix CV= e(scpccvs)
matrix S= e(scpcstats)

matrix list CV
matrix list S
