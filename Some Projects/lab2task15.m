alphas = 1:1:10;
colors = ["black", "red", "blue", "green", "yellow", "black", "red", "blue", "green", "yellow"];
edges = ["none", "blue", "green", "red", "black", "none", "blue", "green", "red", "black"];

drawManyBalls(alphas, colors, edges);

function [] = drawManyBalls(alphas, colors, edges)
    sz = size(alphas, 2);    
    
    f = @(X, Y, Z, i) (abs(X).^alphas(i)+abs(Y).^alphas(i)+abs(Z).^alphas(i)).^(1./alphas(i));
    a = -1;
    b = 1;
    h = 0.1;
    [X,Y,Z] = meshgrid(a:h:b, a:h:b, a:h:b);
    FAlpha = 1:-1/sz:0.1; 
    for i = 1:sz
        F = f(X,Y,Z,i);
        g = isosurface(X,Y,Z,F,1);
        if(size(g.vertices) == 0)
            error('Empty')
        end
        p = patch(g);
        isonormals(X,Y,Z,F,p);
        p.FaceColor = colors(i);
        p.EdgeColor = edges(i);
        p.FaceAlpha = FAlpha(i);
    end
    daspect([1 1 1]);
    view(3);
    axis tight;
    camlight;
    lighting gouraud;
end






