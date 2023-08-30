function [W,delta,n_dropped] = get_PC(r,s,q_in,o,drop_flag)
    
    [S,delta,n_dropped] = get_bspline_mat(s,q_in,o,drop_flag);
    [~,q_max] = size(S);
    q =  min([r,q_max]);
    SS = S'*S;
    [V,D] = eig(SS,'vector');
    [~,index] = sort(D,'descend');
    V_sorted = V(:,index(1:q));
    W = S*V_sorted;

end