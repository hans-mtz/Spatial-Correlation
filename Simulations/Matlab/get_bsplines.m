function [bspline_matrix, delta, n_dropped] = get_bsplines(s,nodes, order)
    [T,n_s] = size(s);
    if n_s > 1
        old = ones(T,1);       
        delta_v = [];
        for i=1:n_s
            new = [];
            [bspline_temp, delta_si, ~] = get_bspline_mat(s(:,i),nodes,order,0);
            for j=1:size(bspline_temp,2)
                new = [new bspline_temp(:,j).*old];
            end
            old = new;
            delta_v = [delta_v delta_si];
        end
        bspline_matrix = new;
        % Dropping columns with no obs to avoid colinearity
        select=sum(bspline_matrix,1)<=1;
        n_dropped=sum(select);
        bspline_matrix(:,select)=[];
        delta = max(delta_v); %Delta should be the same for all
    else    
        [bspline_matrix, delta, n_dropped]=get_bspline_mat(s,nodes,order,1);
    end
end