clear; clc;
f = @(x) x^2 * (2 + 3*x*sin(x) + 4*x^2*(sin(x))^2) - 1;
x0 = [0.2, 0.6];
r = fzero(f, x0);

J = @(x) -4*x - 9*x^2*sin(x) - 3*x^3*cos(x) - 16*x^3*(sin(x))^2 ...
    - 4*x^4*sin(2*x);

J1 = J(r);
J2 = J(-r);

A1 = [0, 1; J1, 1];
A2 = [0, 1; J2, 1];
D1 = eig(A1);
D2 = eig(A2);