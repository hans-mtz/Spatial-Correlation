function epsilon = get_circle_location_noise(s, percent)
    [T, N_s] = size(s);
    i=1;
    radius = percent*1;
    epsilon = NaN(T,N_s);
    while i<T+1
        e = -radius+2*radius*rand(1,N_s);
        if norm(e,2)<=radius
            epsilon(i,:)=e;
            i=i+1;
        end
    end
end