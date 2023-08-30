job17=batch('optimal_spline_sims')
myCluster = parcluster('local');
job14 = findJob(myCluster, 'ID', 14);
diary(job17)