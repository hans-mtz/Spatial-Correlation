# Replication code for Monte Carlo Simulation

Authors: Jianfei Cao, Christian Hansen, Damian Kozbur, Lucciano Villacorta

Date: July 16, 2021

These are the replication files for the Monte Carlo simulation in "Inference for Dependent Data with Learned Clusters" by Jianfei Cao, Christian Hansen, Damian Kozbur, and Lucciano Villacorta.
All codes are written in Matlab and all simulation results are obtained under Matlab R2017b (9.3.0.713579). 

## Structure
There are 8 folders, each of which generates results for one setting:

| Folder Name | Model | Spatial Model | Number of locations | Approximate computation time (s) |
| --- | --- | --- | --- | --- |
| sim1\_ols\_small_baseline | OLS	 | BASELINE | 205 | 5756 |
| sim2\_ols\_small_sar | OLS | SAR | 205 | 6342 |
| sim3\_iv\_small_baseline | IV | BASELINE | 205 | 9485 | 
| sim4\_iv\_small_sar | IV | SAR | 205 | 9245 | 
| sim5\_ols\_large_baseline | OLS | BASELINE | 820 | 33332 |
| sim6\_ols\_large_sar | OLS | SAR | 820 | 27526 | 
| sim7\_iv\_large_baseline | IV | BASELINE | 820 | 36215 |
| sim8\_iv\_large_sar | IV | SAR | 820 | 37318 |

All simulation results in the main text are included in sim1-4. 
Extra simulation results with a larger set of locations (reported in the supplementary materialS) are included in sim5-8. 

For each setting, there are four subfolders, /code/, /data/, /output/, and /temp/. 
To generate the results for a specific setting, run /code/main.m under the corresponding folder. The results (5 csv tables and 1 pdf figure) are reported under /output/. 
The location data from the empirical application is stored under /data/ and temporary files are under /temp/. 

The run_all.m file generates output for all eight settings, although we suggest run each setting separately.  

## Replicating tables and figures

To replicate tables and figures in the paper, check the needed output files in the following table:

| Table | Simulation No. | Output file under /output/ |
| --- | --- | --- |
| Table 3: Simulation Results | sim1-4 | table\_simulation_results.csv |
| Table 4: Distribution of $\widehat{k}$ | sim1-4 | table\_kHat_distribution.csv |
| Table 5: Distribution of $\widehat{\alpha}$ | sim1-4 | table\_alphaHat\_distribution.csv |
| Table 6: Summary: OLS - BASELINE ($N_{\rm{pan}}=205$) | sim1 | table\_summary\_complete.pdf |
| Table 7: Summary: OLS - SAR ($N_{\rm{pan}}=205$) | sim2 | table\_summary\_complete.pdf |
| Table 8: Summary: IV - BASELINE ($N_{\rm{pan}}=205$) | sim3 | table\_summary\_complete.pdf |
| Table 9: Summary: IV - SAR ($N_{\rm{pan}}=205$) | sim4 | table\_summary\_complete.pdf |
| Table 10: Clustering: OLS - BASELINE ($N_{\rm{pan}}=205$) | sim1 | table\_clustering\_complete.pdf |
| Table 11: Clustering: OLS - SAR ($N_{\rm{pan}}=205$) | sim2 | table\_clustering\_complete.pdf |
| Table 12: Clustering: IV - BASELINE ($N_{\rm{pan}}=205$) | sim3 | table\_clustering\_complete.pdf |
| Table 13: Clustering: IV - SAR ($N_{\rm{pan}}=205$) | sim4 | table\_clustering\_complete.pdf |
| Table 14: Summary: OLS - BASELINE ($N_{\rm{pan}}=820$) | sim5 | table\_summary\_complete.pdf |
| Table 15: Summary: OLS - SAR ($N_{\rm{pan}}=820$) | sim6 | table\_summary\_complete.pdf |
| Table 16: Summary: IV - BASELINE ($N_{\rm{pan}}=820$) | sim7 | table\_summary\_complete.pdf |
| Table 17: Summary: IV - SAR ($N_{\rm{pan}}=820$) | sim8 | table\_summary\_complete.pdf |
| Table 18: Clustering: OLS - BASELINE ($N_{\rm{pan}}= 820 $) | sim5 | table\_clustering\_complete.pdf |
| Table 19: Clustering: OLS - SAR ($N_{\rm{pan}}= 820 $) | sim6 | table\_clustering\_complete.pdf |
| Table 20: Clustering: IV - BASELINE ($N_{\rm{pan}}= 820 $) | sim7 | table\_clustering\_complete.pdf |
| Table 21: Clustering: IV - SAR ($N_{\rm{pan}}= 820 $) | sim8 | table\_clustering\_complete.pdf |

| Figure | Simulation No. | Output file under /output/ |
| --- | --- | --- |
| Figure 2: OLS power curves | sim1-2 | power\_curve\_ols\_baseline\_small.pdf, power\_curve\_ols\_sar\_small.pdf | 
| Figure 3: IV power curves | sim3-4 | power\_curve\_iv\_baseline\_small.pdf, power\_curve\_iv\_sar_small.pdf |
| Figure 4: OLS power curves ($N_{\rm{pan}}=820$) | sim5-6 | power\_curve\_ols\_baseline\_large.pdf, power\_curve\_ols\_sar\_small.pdf | 
| Figure 5: IV power curves ($N_{\rm{pan}}=820$) | sim7-8 | power\_curve\_iv\_baseline\_large.pdf, power\_curve\_iv\_sar_small.pdf |




