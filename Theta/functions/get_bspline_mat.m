function [S, delta, n_dropped] = get_bspline_mat(x,n_nodes,order,drop_flag)
    % O = order;
    n = length(x);
    a = 0; %min(x);
    b = 1; %max(x);
    % if order <=2
    %     delta = (b-a)/(n_nodes-1);
    % else
    %     delta = (b-a)*3/(n_nodes-2)*2;
    % end

    if order==1
        delta = (b-a)/(n_nodes+1);
        t = linspace(a,b,n_nodes+1);
    elseif order ==2
        delta = (b-a)/(n_nodes-1);
        t = linspace(a-delta,b+delta,n_nodes+2);
    elseif order ==3
        delta = (b-a)/(n_nodes-2);
        t = linspace(a-delta*(3),b+delta*(3),n_nodes+3);
    end
    
    % fprintf('t %d\n',t)
    % fprintf('delta %d\n',delta)
    T = length(t);

    B = bspline1(t,order);

    bspline_mat = [];
    
    for i=1:n
        b_row = [];
        for j=order %1:O Only get the Splines for the order we want
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
    
    % if K>1
    %     S(:,1:T-1) = [];
    % else
    %     S(:,[1 T-1]) = [];
    % end
%     S(:,[1,T-1])=[]; % Removing columns outside the range [1,T-1] bc they are only zeros
    % Removing one dummy column, this is the reference dummy to avoid
    % multicolinearity
    if drop_flag==1
%         if K ==1
%             S(:,2) = [];
%         end
        select=sum(S,1)<=1;
        n_dropped=sum(select);
        S(:,select)=[];
    else
        n_dropped=0;
    end
end