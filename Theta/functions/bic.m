function value = bic(residuals,X,method)
    [n, k] = size(X);
    % SSR=residuals'*residuals;
    sigma_2 = var(residuals);
    switch method
        case 'hansen'
            value = n*log(2*pi*sigma_2)+n+(k+1)*log(n);      
        otherwise
            value = n*log(sigma_2) + k*log(n);
    end
end