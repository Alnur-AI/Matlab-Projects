function res = Regression(n,x,y,k)
    %n размеры x,y. k-степень многочлена
    A = zeros(k+1,k+1);
    s = zeros(1,2*k+1);
    s(1) = n;
    for i = 2:2*k
        s(i)=sum((x.^(i-1)));
    end
    for i=1:k+1
        for j=1:k+1
            A(i,j)=s(i+j-1);
        end
    end
    z = zeros(k+1,1);
    z(1) = sum(y);
    for (i=2:k+1)
        z(i) = sum(x.^(i-1).*y);
    end
    res = linsolve(A,z);
    
end