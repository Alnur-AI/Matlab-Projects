h = 0.5;
a = -2;
b = 2;
alpha = 2;
f = @(x, y, z) abs(x) .^ alpha + abs(y) .^ alpha + abs(z) .^ alpha;

params.f = f;
params.a = a;
params.b = b;
params.h = h;
params.isovalue = 100;
params.FaceColor = 'red';
params.EdgeColor = 'blue';
%params.color = 'green';
params.FaceAlpha = 1;
params.LineStyle = '-';
params.Marker = 'o';
params.MarkerFaceColor = 'green';
params.MarkerSize = 5; 


drawBall(alpha, params);
%%
function [] = drawBall(alpha, params)
    [X, Y, Z] = meshgrid(params.a:params.h:params.b, params.a:params.h:params.b, params.a:params.h:params.b);
    F = params.f(X,Y,Z); 
    g = isosurface(X,Y,Z,F,params.isovalue)    
    p = patch(g);
    if(size(g.vertices) == 0)
        error('Empty')      
    end
    isonormals(X,Y,Z,F,p);
    view(3);
    p.FaceColor = params.FaceColor;
    p.EdgeColor = params.EdgeColor;
    p.FaceAlpha = params.FaceAlpha;
    p.LineStyle = params.LineStyle;
    p.Marker = params.Marker;
    p.MarkerFaceColor = params.MarkerFaceColor;
    p.MarkerSize = params.MarkerSize;
    daspect([1 1 1]);
    axis tight;
    camlight;
    lighting gouraud;
end


