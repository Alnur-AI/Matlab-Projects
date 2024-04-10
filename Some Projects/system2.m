function dxdt = system2(t, x)
    f = @(a) sign(x(4)) - 2*a^2 - 3*a^3*sin(a) - 4*a^4*(sin(a))^2;
    g = @(a) - 4*a - 9*a^2*sin(a) - 3*a^3*cos(a) - 16*a^3*(sin(a))^2 ...
             - 4*a^4*sin(2*a);
    dxdt = [x(2); f(x(1)) + x(2); x(4) * g(x(1)); x(3) + x(4)];
end