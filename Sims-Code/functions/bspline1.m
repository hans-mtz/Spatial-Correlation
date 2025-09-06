function B =bspline1(t,K)

T=length(t);
B = cell(T-1,K);

for i=1:T-1
    B{i,1} = @(x) abs(x>=t(i) & x < t(i+1));
end
for k=2:K
    for i=1:T-1
        if i+k <= T
            B{i,k} = @(x) (x-t(i))/(t(i+k-1)-t(i)).*B{i,k-1}(x) ...
                        + (t(i+k)-x)/(t(i+k)-t(i+1)).*B{i+1,k-1}(x);
        end
    end
end

end





















