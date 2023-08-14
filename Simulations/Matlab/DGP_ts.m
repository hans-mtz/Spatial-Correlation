function [y,X,u] = DGP_ts(beta,coef,T, model)
    %rho_bar 0.03 -> AR1 0.79
    burn = 100;
    k = length(beta);
    u_temp = zeros(T+burn,1);
    x_temp = zeros(T+burn,1);
    epsilon = randn(T+burn,2);
    for i=1:T+burn
        u_temp(i+1) = u_temp(i)*coef + epsilon(i,1);
        x_temp(i+1) = x_temp(i)*coef + epsilon(i,2);
    end
    u = u_temp(burn+1:end-1);
    % size(u)
    switch model
        case 1
            X = [ones(T,1) randn(T,k-1)];
        otherwise
            X = [ones(T,1) x_temp(burn+1:end-1)];
    end
    y = X*beta + u;

end