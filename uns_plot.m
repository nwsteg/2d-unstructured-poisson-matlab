function uns_plot( cells )

addpath('./meshlib');

num=length(cells);
X=zeros(num,1);
Y=zeros(num,1);
for c = cells
    X(c.id)=c.cen.x;
    Y(c.id)=c.cen.y;
end

tri=delaunay(X,Y);
trisurf(tri,X,Y,ones(num_cells,1));
axis vis3d;

end

