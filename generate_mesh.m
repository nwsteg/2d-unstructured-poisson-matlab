function [ vert,etri,tria,tnum,node,edge ] = generate_mesh( hfun,plot )

addpath('./meshlib');

% Define a 1x1 square
node = [                % list of xy "node" coordinates
    0, 0                % outer square
    1, 0
    1, 1
    0, 1 ] ;

edge = [                % list of "edges" between nodes
    1, 2                % outer square 
    2, 3
    3, 4
    4, 1 ] ;

% call mesh generator
% vert -- the vertices in the mesh
% etri -- 
[vert,etri, ...
tria,tnum] = refine2(node,edge) ;

[vert,etri, ...
tria,tnum] = refine2(node,edge,[],[],hfun) ;

if(plot==1)
    figure;
    patch('faces',tria(:,1:3),'vertices',vert, ...
        'facecolor','w', ...
        'edgecolor',[.2,.2,.2]) ;
    hold on; axis image off;
    patch('faces',edge(:,1:2),'vertices',node, ...
        'facecolor','w', ...
        'edgecolor',[.1,.1,.1], ...
        'linewidth',1.5) ;
end

end

