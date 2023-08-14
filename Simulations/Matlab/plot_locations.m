%% Plotting locations and distances %%

% setting seed
rng(333)
T = 2500;
B = 20;
% Uniformly distributed locations

s_u = rand(T,1); % vector of locations

% Time series

s_ts = (1:T)'./T;

D_u = getdistmat(s_u,0);
D_ts = getdistmat(s_ts,0);

%% Plotting locations %%

h=histogram(s_u,'NumBins',B);
h_val=h.Values;
h_bin=h.BinEdges;

h_ts=histogram(s_ts,'NumBins',B);
h_ts_val=h_ts.Values;
% h_ts_bin=h_ts.BinEdges;

bar(h_bin(2:end),[h_val' h_ts_val'])
exportgraphics(gcf,strcat('../Products/hist_loc_',num2str(T),'.png'))

%% Plotting distances %%

h_dist=histogram(D_u,'NumBins',B);
h_dist_val=h_dist.Values;
h_dist_bin=h_dist.BinEdges;

h_dist_ts=histogram(D_ts,'NumBins',B);
h_dist_ts_val=h_dist_ts.Values;
% h_dist_ts_bin=h_dist_ts.BinEdges;

bar(h_dist_bin(2:end),[h_dist_val' h_dist_ts_val'])
exportgraphics(gcf,strcat('../Products/hist_dist_',num2str(T),'.png'))

