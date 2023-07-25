function [y,X,D_mat] = DGP(beta,T,k,s,rho_bar,model)

% Generate variables
% Z = randn(T,k); % matrix of exogenous variables
D_mat = getdistmat(s,false); % matrix of distances
c = getcbar(rho_bar, D_mat); % calculating the c_min given rho_bar
Sigma = exp(-c*D_mat); % var-cov matrix 
% rho_estimated = mean(reshape(Sigma,[],1)); % estimated c_min corresponds to 
% fprintf('Estimated c_min corresponds to a rho_bar estimated of %5.4g \n', rho_estimated);
% mu = zeros(T,1); % mean zero for error terms
% u = mvnrnd(mu,Sigma)'; % sampling from G_exp(c_min)

switch model
    case 2
        mu_x = zeros(k-1,T);
        X_mat = mvnrnd(mu_x,Sigma);
        X = [ones(T,1) X_mat'];
    otherwise
        X = [ones(T,1) randn(T,k-1)]; %iid
end

switch rho_bar
    case 0.0
        u = randn(T,1);
    otherwise
        mu = zeros(T,1); % mean zero for error terms
        u = mvnrnd(mu,Sigma)'; % sampling from G_exp(c_min)
end

y = X*beta + u;

% dt.y = y;
% dt.X = X;
% dt.D_mat = D_mat;

end