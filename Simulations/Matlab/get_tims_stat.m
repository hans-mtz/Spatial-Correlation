function var_stat = get_tims_stat(u)
    % T = size(u,1);
    % mu = sum(u)/T;
    % VC = (u-mu)*(u-mu)';
    VC = u*u';

    % var_stat = (sum(VC-diag(diag(VC)),'all')/(T*(T-2)-1))/(sum(diag(VC))/(T*(T-1)-1));
    var_stat = sum(VC-diag(diag(VC)),'all')/sum(diag(VC));
end