function I = my_simpson(x,y)
    
    n = sum( size(x) ) - 2;
    a = x(1,1);
    b = x(1,end);
    h = (b-a)/n;

    if mod(n+1,2) == 0
        %error('n is odd!');
        n = n+1;
        h = (b-a)/n;
        
        x = a:h:b; %n+1 elements
        y = dSi(x);%n+1 elements
        
        
        y1 = y(1,1:2:end-2);
        y2 = y(1,2:2:end-1);
        y3 = y(1,3:2:end);

        I = sum( (y1(:) + 4*y2(:) + y3(:)) );

        I = I*(h/3);
        return;
    end
    
    
    
    I = 0;
    
    y1 = y(1,1:2:end-2);
    y2 = y(1,2:2:end-1);
    y3 = y(1,3:2:end);
    
    I = sum( (y1(:) + 4*y2(:) + y3(:)) );
    
    I = I*(h/3);
    
end
%{
 I = 0;
        n = n+1;
        h = (b-a)/n;
        
        x = a:h:b; %n+1 elements
        y = dSi(x);%n+1 elements
        
        
        y1 = y(1,1:2:end-2);
        y2 = y(1,2:2:end-1);
        y3 = y(1,3:2:end);

        I = sum( (y1(:) + 4*y2(:) + y3(:)) );

        I = I*(h/3);
        return;
%}