function [S] = get_bspline_mat(x,n_nodes,order)
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
    S(:,1:T-1) = []; 
%     S(:,[1,2,T-1])=[]; % Removing columns outside the range [1,T-1] bc they are only zeros
    % Removing one dummy column, this is the reference dummy to avoid
    % multicolinearity
%     S = bspline_mat(:,2:(size(bspline_mat,2)-1));
end