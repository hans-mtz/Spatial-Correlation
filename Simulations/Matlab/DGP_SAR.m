function [y,X,u] = DGP_SAR(beta,coef,D_mat,L)
    T = size(D_mat,1);
    k = length(beta);
    W = D_mat < L;
    W = W - diag(diag(W));

    epsilon = randn(T,k);

    u = (eye(T)-coef*W)\epsilon(:,1);
    X = [ones(T,1) (eye(T)-coef*W)\epsilon(:,2:k)];

    y=X*beta+u;

end