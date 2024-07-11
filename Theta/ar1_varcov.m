if batchStartupOptionUsed
    addpath(genpath('./functions'))
    % addpath(genpath('./21200057'))
end

rng(549)


%% Setting up parameters %%%%

n_reps = 200;
T = 500;
n_locations = 2;
theta = sqrt(2)/10;


%% Generate locations (fixed locations)

s = rand(T,n_locations);
delta=[0.94:0.01:0.97];

A=createMatrix(T-1);
b=[];
for d=delta
    b=[b sum(d.^A,'all')/(T*T)];
end
% c=[]
% for i=1:200
%     [y, X, h] = DGP(theta,s,0.75);
%     c=[c sum(X(:,2)*X(:,2)',"all")/(n*n)];
% end
% mean(c)

[y, X, h] = DGP(theta,s,0.75);

sum(exp(-h/theta),"all")/(T*T)

array2table([delta' b' ],"VariableNames",["delta" "sum"])