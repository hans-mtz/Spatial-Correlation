function [y,X,D_mat] = DGP(theta,locations,rho,intercept,D_mat)
    arguments
        theta (1,1) double
        locations (:,:) double
        rho (1,:) double
        intercept (1,1) logical = true
        D_mat (:,:) double = getdistmat(locations,false)
    end

    T=length(locations);

    % D_mat = getdistmat(locations,false);

    Sigma = (1-rho).*eye(T) +rho.*exp(-D_mat/theta);
    mu = zeros(2,T);
    Z = mvnrnd(mu,Sigma);
    if (intercept)
        X = [ones(T,1) Z(1,:)'];
    else 
        X = Z(1,:)';
    end
    y = Z(2,:)';

    % % Spatial correlation
    % % Sigma = exp(-D_mat/theta);
    % mu = zeros(2,T);
    % X_spatial = mvnrnd(mu,Sigma);

    % % White noise
    % X_noise = normrnd(0,1,[2,T]);
    

    % % X= [ones(T,1) (1-rho).*X_noise(1,:)'+rho.*X_spatial(1,:)']; % With Intercept
    % X= (1-rho).*X_noise(1,:)'+rho.*X_spatial(1,:)'; %No Intercept

    % % Y is an independent draw of the same GDP

    % y = (1-rho).*X_noise(2,:)'+rho.*X_spatial(2,:)';

end