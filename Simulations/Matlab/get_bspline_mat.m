function [S, delta, n_dropped] = get_bspline_mat(x,n_nodes,order,drop_flag)
    K = order;
    n = length(x);
    a = 0; %min(x);
    b = 1; %max(x);
    delta = (b-a)/(n_nodes-1);
    t = linspace(a-delta,b+delta,n_nodes+2);
    T = length(t);

    B = bspline1(t,K);

    
    bspline_mat = [];
    for i=1:n
        b_row = [];
        for j=1:K
            for k=1:T-1
                if j+k <= T
%                     display([i,k,j]);
                    b_row = [b_row B{k,j}(x(i))];
                end
            end
        end
%         display(b_row);
        bspline_mat = [bspline_mat; b_row];
    end
    
    S = bspline_mat;
    
    if K>1
        S(:,1:T-1) = [];
    else
        S(:,[1 2 T-1]) = [];
    end
%     S(:,[1,T-1])=[]; % Removing columns outside the range [1,T-1] bc they are only zeros
    % Removing one dummy column, this is the reference dummy to avoid
    % multicolinearity
    if drop_flag==1
        select=sum(S,1)<=1;
        n_dropped=sum(select);
        S(:,select)=[];
    else
        n_dropped=0;
    end
end