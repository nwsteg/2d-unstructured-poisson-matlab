% plot the result
function plot_uns(X,Y,u)
    tri=delaunay(X,Y);
    [r,c]=size(tri);
    h=trisurf(tri,X,Y,u);
%     shading interp;
%     axis vis3d;
end