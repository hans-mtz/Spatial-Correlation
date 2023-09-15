function nn_corr = get_nn_corr(u,D,L)
    T = size(u,1);
    if length(L)>1
        nn_corr = NaN(length(L),1);
    else
        nn_corr = NaN;
    end
    
    mu = mean(u*u'-diag(diag(u*u')),'all');
    k=1;
    for l=L
        rho_i=NaN(T,1);
    for i=1:T
        rho=[];%NaN(T,T);
        if length(L)>1
            index=find(D(i,:)<=l);
        else
            temp = D(i,1:end ~= i);
            val = min(temp);
            index = find(D(i,:)<=val);
        end

        for j=index
            if i~=j
                rho = [rho (u(i)*u(j)-mu)];

            end
        end
    %     rho_i(i) = (mean(rho))/(sigma);
    %     rho_i(i) = (mean(rho)-(mu*mu))/norm(u_hat(index));
        rho_i(i) = mean(rho);
    end
    nn_corr(k) = sum(rho_i)./((T-1));
    k=k+1;
    end
end