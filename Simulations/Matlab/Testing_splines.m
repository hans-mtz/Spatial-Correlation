% Testing BSplines %
K = 2;
x = 0.0:0.001:1.0;
n = length(x);
% t = linspace(0.0,1.0,4);
t = linspace(-0.5,1.5,6); %B-Spline space should be
% greater to get the half trianges at the edges
T = length(t);

B = bspline1(t,K);

B_mat = NaN(n,(T-1),K);
% J = size(B,2);
for i=1:n
    for k=1:T-1
        for j=1:K
            if j+k <= T
                B_mat(i,k,j) = B{k,j}(x(i));
            end
        end
    end
end

plot(x,B_mat(:,:,1)); %B Spline order 1 -> Step function
plot(x,B_mat(:,:,2)); %B Spline order 2 -> Triangles
plot(x,[B_mat(:,1:2,2) B_mat(:,:,1)]);

% get_bspline_mat function works
x = 0.0:0.001:1.0;
A=get_bspline_mat(x,4,1); %
plot(x,A);

x = 0.0:0.001:1.0;
S=get_bspline_mat(s,125,2); %
sum(S,1)
sum(S,2)
plot(x,A);

x = 0.0:0.001:1.0;
A=get_bspline_mat(x,120,1); %

plot(x,A);

x = 0.0:0.001:1.0;


for j=[4 5 6 7]
    A=get_bspline_mat(x,j,2); %
    plot(x,A);
    exportgraphics(gcf,strcat('../Products/bs_plot_',num2str(j),'.png'))
end

S=get_bspline_mat(s,4,1);
S
ols(y,[X(:,2:end) S],[X(:,2:end) S], 'ldl');
ols(y,[X(:,2:end) S],[X(:,2:end) S], 'chol');
ols(y,[X S],[X S], 'ldl');
ols(y,[X S],[X S], 'chol');
ols(y,
