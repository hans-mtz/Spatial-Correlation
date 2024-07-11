function [y,X,D_mat] = DGP(theta,locations,rho)
    T=length(locations);

    D_mat = getdistmat(locations,false);

    % Spatial correlation
    Sigma = exp(-D_mat/theta);
    mu = zeros(2,T);
    X_spatial = mvnrnd(mu,Sigma);

    % White noise
    X_noise = normrnd(0,1,[2,T]);

    X= [ones(T,1) (1-rho).*X_noise(1,:)'+rho.*X_spatial(1,:)'];

    % Y is an independent draw of the same GDP

    y = (1-rho).*X_noise(2,:)'+rho.*X_spatial(2,:)';

end