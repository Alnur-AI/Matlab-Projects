function C = my_multiply(A,B)
    As = size(A);
    Bs = size(B);
    
    C = zeros(As(1,1), Bs(1,2));
    
    for i = 1:As(1,1)
        for j = 1:Bs(1,2)
            for k = 1:Bs(1,1)
                C(i,j) = C(i,j) + A(i,k)*B(k,j);
            end
        end
    end
    
end

