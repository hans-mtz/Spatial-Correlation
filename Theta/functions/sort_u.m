function [u_t, u_t_1] = sort_u(u,s)
    T = length(s);
    [~, sorting_index] = sort(s); % Sorting by locations
    u_sorted = u(sorting_index);
    u_t = u_sorted(1:end-1);
    u_t_1 =  u_sorted(2:end);
    u_t_1 = [ones(T-1,1) u_t_1];
end