function y = c_derivative(x,h)
    y = (f13(x+h) - f13(x-h))./(2*h);
end
