function [X, Y, Sx, Sy] = reachset(T, whole, showTrace)
    [X, Y, Sx, Sy] = createSide(1, 1e-8, 1e-8, T, showTrace);
    if (whole == 1)
        [X2, Y2] = createSide(-1, 1e-8, 1e-8, T, showTrace);
        X = cat(1, X, X2);
        Y = cat(1, Y, Y2);
    end
end

function [X, Y, Sx, Sy] = createSide(u, absTol, relTol, T, showTrace)
    options = odeset('Events', @zeroX, 'RelTol',relTol,'AbsTol',absTol);

    t0 = 0;
    x0 = [0, 0];
    f1 = @(t, x) system1(t, x, u);
    [t1, x1] = ode45(f1, [t0 T], x0, options);
    
    tau1 = t1(end);
    N = size(t1, 1);

    options = odeset('RelTol',relTol,'AbsTol',absTol);
    Sx = x1(:, 1);
    Sy = x1(:, 2);
%save loverest line
    i = 2;
    tau = t1(i);
    x01 = [x1(i, 1); x1(i, 2); sign(x1(i, 2)); (-u)*1e-3];
    [t, x] = ode45(@system2, [tau T], x01, options);
    X = x(:, 1);
    Y = x(:, 2);
%save ends of other lines
    for i = 3 : N-2
        tau = t1(i);
        x01 = [x1(i, 1); x1(i, 2); sign(x1(i, 2)); (-u)*1e-4];
        [t, x] = ode45(@system2, [tau T], x01, options);
        if (showTrace == 1)
            plot(x(:, 1), x(:, 2), '-r');%debug plot
        end
        X(size(X, 1) + 1) = x(end, 1);
        Y(size(Y, 1) + 1) = x(end, 2);
    end
    i = i + 1;
    tau = t1(i);
    x01 = [x1(i, 1); x1(i, 2); sign(x1(i, 2)); (-u)*1e-4];
    [t, x] = ode45(@system2, [tau T], x01, options);
    
    if (showTrace == 1)
        plot(x(:, 1), x(:, 2), '-r');%debug plot
    end
    
    X = cat(1, X, flip(x(:, 1)));
    Y = cat(1, Y, flip(x(:, 2)));
    
    X = cat(1, X, flip(x1(:, 1)));
    Y = cat(1, Y, flip(x1(:, 2)));

end
















