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

%% AR(1)
delta=[0.90:0.01:0.99];
T=500
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

tbl_500=array2table([delta' b' ],"VariableNames",["delta" "sum"])

%% AR(1)
% delta=[0.90:0.01:0.99];
T=250;
s = rand(T,n_locations);
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

tbl_250=array2table([delta' b' ],"VariableNames",["delta" "sum"])

%% Spatial

[y, X, h] = DGP(theta,s,0.75);

sum(exp(-h/theta),"all")/(T*T)

