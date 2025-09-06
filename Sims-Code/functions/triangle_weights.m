function K = triangle_weights(s,l)
    [T, d] = size(s);

    % Get distances per locations dimention
    D_mat_dim_pairs= NaN(T,T,d);
    for i=1:d
        % tmp = repmat(s(:,i),1,2);
        % D_mat_dim_pairs(:,:,i) = squareform(pdist(tmp));
        % D_mat_dim_pairs(:,:,i) = getdistmat(tmp,false);
        D_mat_dim_pairs(:,:,i) = getdistmat(s(:,i),false);
    end 

    % Get weight per dimention of location

    K_temp = zeros(T,T,d);

    for i=1:T
        for dim=1:d
            for j=find(D_mat_dim_pairs(i,:,dim)<l)
                K_temp(i,j,dim) = triangle(D_mat_dim_pairs(i,j,dim),l);
            end
        end
    end

    K = ones(T,T);
    for dim=1:d
        K(:,:) = K.*K_temp(:,:,dim);
    end

end

