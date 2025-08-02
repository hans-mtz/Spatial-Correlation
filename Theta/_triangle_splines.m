
%% Setting directory and adding functions to path
% cd Theta
addpath(genpath('./functions'))

%% B Splines
t= 0:0.005:1;
S = get_bsplines(t',8,2);

%% Plot
scatter(t,S,[], "filled")
title('8 Triangle Spline, Illustration')
gray_scale = [0 0 0; 
              0.1 0.1 0.1; 
              0.2 0.2 0.2; 
              0.3 0.3 0.3; 
              0.4 0.4 0.4; 
              0.5 0.5 0.5; 
              0.6 0.6 0.6; 
              0.7 0.7 0.7];
colororder(gray_scale)
exportgraphics(gcf,strcat(['figures/filled_circles_splines_plot.pdf']))
%% Plot
scatter(t,S)
title('8 Triangle Spline, Illustration')
gray_scale = [0 0 0; 
              0.1 0.1 0.1; 
              0.2 0.2 0.2; 
              0.3 0.3 0.3; 
              0.4 0.4 0.4; 
              0.5 0.5 0.5; 
              0.6 0.6 0.6; 
              0.7 0.7 0.7];
colororder(gray_scale)
exportgraphics(gcf,strcat(['figures/empty_circles_splines_plot.pdf']))
%% B Splines
t= 0:0.001:1;
S = get_bsplines(t',8,2);


%% lines
plot(t,S)
title('8 Triangle Spline, Illustration')
exportgraphics(gcf,strcat(['figures/lines_splines_plot.pdf']))