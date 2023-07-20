clear all

import delimited "/Volumes/SSD Hans 1/Github/Spatial Correlation/Simulations/Products/example.csv"

// Regression
reg y x1, nocon robust

// Spatial correlation comand
scpc, cvs
