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