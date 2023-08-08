function [y,X,D_mat] = DGP(beta,s,rho_bar,model)
% DGP - M1 X~{1 iid}: M2 X~{1 G_exp}; M3 M2 & Y demeaned;

T=length(s);
k=length(beta);
% Generate variables

D_mat = getdistmat(s,false); % matrix of distances
c = getcbar(rho_bar, D_mat); % calculating the c_min given rho_bar
Sigma = exp(-c*D_mat); % var-cov matrix


% Correlated regressors

mu_x = zeros(k-1,T);
X_mat = mvnrnd(mu_x,Sigma);

% error term
switch rho_bar
    case 0.0 % iid
        u = randn(T,1);
    otherwise % spatially correlated 
        mu = zeros(T,1); % mean zero for error terms
        u = mvnrnd(mu,Sigma)'; % sampling from G_exp(c_min)
end

% Models
switch model
    case 3
        X = [ones(T,1) X_mat'];
        y_0 = X*beta + u;
        y = y_0 - mean(y_0);
    case 2 % Spatially correlated
        X = [ones(T,1) X_mat'];
        y = X*beta + u;
    otherwise %iid
        X = [ones(T,1) randn(T,k-1)];
        y = X*beta + u;
end

% y = X*beta + u;


end