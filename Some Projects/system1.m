function dxdt = system1(~, x, u)
    f = @(a) u - 2*a^2 - 3*a^3*sin(a) - 4*a^4*(sin(a))^2;
    dxdt = [x(2); f(x(1)) + x(2)];
end