function I = my_rectangles(x,y)
    
    n = sum( size(x) ) - 2;
    a = x(1,1);
    b = x(1,end);
    h = (b-a)/n;
    
    
    I = 0;
    
    x1 = x(1,2:end);
    x2 = x(1,1:end-1);
    y1 = y(1,2:end);
    I = sum( y1(:).*(x1(:) - x2(:)) );
end
