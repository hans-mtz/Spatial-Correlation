%% Loading %%
clear;
excercise = 'dd_2x';%'grid';%

load(['bic_spline_' excercise '.mat'])

%% Plotting %%
fprintf('Plotting\n');
nbins=50;
histogram(splines_array(:,1),nbins,'BinLimits',[10,90]);
exportgraphics(gcf,strcat(['../Products/hist_n_splines_o1_bic_' excercise '.png']))

histogram(splines_array(:,2),nbins,'BinLimits',[10,90]);
exportgraphics(gcf,strcat(['../Products/hist_n_splines_o2_bic_' excercise '.png']))

%% Double bar graph %%

clear;
excercise = 'dd_2x';%'grid';%

load(['bic_spline_' excercise '.mat'])

splines_dd_2x = splines_array;

clear splines_array;

excercise = 'grid';%

load(['bic_spline_' excercise '.mat'])

splines_grid = splines_array;

%% 

nbins=20;
h_splines_dd_2x = histogram(splines_dd_2x(:,1),nbins,'BinLimits',[10,90]);
h_splines_dd_2x_val = h_splines_dd_2x.Values;
h_splines_dd_2x_BinEdges = h_splines_dd_2x.BinEdges;
% exportgraphics(gcf,strcat(['../Products/hist_n_splines_o1_bic_' excercise '.png']))

h_splines_grid = histogram(splines_grid(:,1),nbins,'BinLimits',[10,90]);
h_splines_grid_val = h_splines_grid.Values;
h_splines_grid_BinEdges = h_splines_grid.BinEdges;
% exportgraphics(gcf,strcat(['../Products/hist_n_splines_o2_bic_' excercise '.png']))
bar(h_splines_dd_2x_BinEdges(1:end-1),[h_splines_dd_2x_val' h_splines_grid_val'])
exportgraphics(gcf,'../Products/hist_n_splines_o1.png')

h_splines_dd_2x = histogram(splines_dd_2x(:,2),nbins,'BinLimits',[10,90]);
h_splines_dd_2x_val = h_splines_dd_2x.Values;
h_splines_dd_2x_BinEdges = h_splines_dd_2x.BinEdges;
% exportgraphics(gcf,strcat(['../Products/hist_n_splines_o1_bic_' excercise '.png']))

h_splines_grid = histogram(splines_grid(:,2),nbins,'BinLimits',[10,90]);
h_splines_grid_val = h_splines_grid.Values;
h_splines_grid_BinEdges = h_splines_grid.BinEdges;
% exportgraphics(gcf,strcat(['../Products/hist_n_splines_o2_bic_' excercise '.png']))
bar(h_splines_dd_2x_BinEdges(1:end-1),[h_splines_dd_2x_val' h_splines_grid_val'])
exportgraphics(gcf,'../Products/hist_n_splines_o2.png')
