function A = createMatrix(n)
    % Inicializar la matriz de ceros
    A = zeros(n+1, n+1);
    
    % Llenar la matriz
    for i = 0:n
        for j = 0:n
            if j >= i
                A(i+1, j+1) = j - i;
            else
                A(i+1, j+1) = i - j;
            end
        end
    end
end