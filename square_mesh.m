% Playing around with MESH2D library
% 6/24/18

% Goal: Figure out how the outputs of refine2 describe
% the resulting mesh.
% Knowing this, I can calculate gradients and define 
% values at the cell centers.

clear all; close all; clc;

addpath('./meshlib');
addpath('./surf_from_scatter');

global cells

origin.x=0; origin.y=0;

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

hfun = +.20; % uniform "target" edge-lengths

[vert,etri, ...
tria,tnum] = refine2(node,edge,[],[],hfun) ;

num_cells=length(tria);

% calculate the triangle centroids
centroid=zeros(num_cells,2);
for i=1:num_cells
    centroid(i,1)=sum(vert(tria(i,:),1))/3;
    centroid(i,2)=sum(vert(tria(i,:),2))/3;
end

% initialize mesh structure
cells = []; % cell
cells.vert=[]; % cell vertices
cells.faces=[];    % cell faces
cells.faces.ss=[]; % cell face start
cells.faces.se=[]; % cell face end
cells.faces.l=0.0; % side length
cells.vol=0.0; % cell volume
cells.cen=[]; % cell centroid


% populate the mesh data structure
for i=1:num_cells
%     fprintf('Triangle %d \n',i);
%     fprintf('Vertices: (%.2f, %.2f), (%.2f, %.2f), (%.2f, %.2f) \n', vert(tria(i,1),1), ... 
%                                vert(tria(i,1),2), ...
%                                vert(tria(i,2),1), ...
%                                vert(tria(i,2),2), ...
%                                vert(tria(i,3),1), ...
%                                vert(tria(i,3),2)  ...
%                                );
%     fprintf('Centroid: (%.2f, %.2f) \n\n',centroid(i,1),centroid(i,2));
%     
    % cell id
    cells(i).id=i;
    
    % cell vertices
    cells(i).vert(1).id=tria(i,1);%cells(i).id*3+1;
    cells(i).vert(1).x=vert(tria(i,1),1);
    cells(i).vert(1).y=vert(tria(i,1),2);
    
    cells(i).vert(2).id=tria(i,2);%cells(i).id*3+2;
    cells(i).vert(2).x=vert(tria(i,2),1);
    cells(i).vert(2).y=vert(tria(i,2),2);
    
    cells(i).vert(3).id=tria(i,3);%cells(i).id*3+3;
    cells(i).vert(3).x=vert(tria(i,3),1);
    cells(i).vert(3).y=vert(tria(i,3),2);
    
    
    % cell faces
    cells(i).faces(1).num=1;
    cells(i).faces(1).id=0;
    cells(i).faces(1).ss=cells(i).vert(1);
    cells(i).faces(1).ss.id=cells(i).vert(1).id;
    cells(i).faces(1).se=cells(i).vert(2);
    cells(i).faces(1).se.id=cells(i).vert(2).id;
    cells(i).faces(1).l=get_length(cells(i).faces(1).ss,...
                                   cells(i).faces(1).se); 
    % face normal
    delx=cells(i).faces(1).ss.x-cells(i).faces(1).se.x;
    dely=cells(i).faces(1).ss.y-cells(i).faces(1).se.y;
    alpha=sqrt(delx*delx+dely*dely);
    cells(i).faces(1).nv.x =  dely/alpha;
    cells(i).faces(1).nv.y = -delx/alpha;
                               
    cells(i).faces(2).num=2;
    cells(i).faces(2).id=0;
    cells(i).faces(2).ss=cells(i).vert(2);
    cells(i).faces(2).ss.id=cells(i).vert(2).id;
    cells(i).faces(2).se=cells(i).vert(3);
    cells(i).faces(2).se.id=cells(i).vert(3).id;
    cells(i).faces(2).l=get_length(cells(i).faces(2).ss,...
                                   cells(i).faces(2).se);
    % face normal
    delx=cells(i).faces(2).ss.x-cells(i).faces(2).se.x;
    dely=cells(i).faces(2).ss.y-cells(i).faces(2).se.y;
    alpha=sqrt(delx*delx+dely*dely);
    cells(i).faces(2).nv.x =  dely/alpha;
    cells(i).faces(2).nv.y = -delx/alpha;
    
    cells(i).faces(3).num=3; 
    cells(i).faces(3).id=0;
    cells(i).faces(3).ss=cells(i).vert(3);
    cells(i).faces(3).ss.id=cells(i).vert(3).id;
    cells(i).faces(3).se=cells(i).vert(1);
    cells(i).faces(3).se.id=cells(i).vert(1).id;
    cells(i).faces(3).l=get_length(cells(i).faces(3).ss,...
                                   cells(i).faces(3).se); 
                               
    % face normal
    delx=cells(i).faces(3).ss.x-cells(i).faces(3).se.x;
    dely=cells(i).faces(3).ss.y-cells(i).faces(3).se.y;
    alpha=sqrt(delx*delx+dely*dely);
    cells(i).faces(3).nv.x =  dely/alpha;
    cells(i).faces(3).nv.y = -delx/alpha;
                           
    % calculate the volume
    cells(i).vol=get_volume(cells(i).vert(1), ...
                            cells(i).vert(2), ...
                            cells(i).vert(3) );
                        
    % calculate the centroid
    cells(i).cen=get_centroid(cells(i).vert(1), ...
                              cells(i).vert(2), ...
                              cells(i).vert(3) );
    
end
% 
% for c = cells
%     fprintf('%d, vol %f, cen (%.2f,%.2f) \n',c.id,c.vol,c.cen.x,c.cen.y);
%     for f = c.faces
%         fprintf('num %d, id %d, ||n|| %.2f \n', f.num, f.id, get_length(f.nv,origin));
%     end
% end

% initialize a field u
u=zeros(num_cells,1);
for c=cells
    x=c.cen.x;
    y=c.cen.y;
    u(c.id)=20*cos(3*pi*x)*sin(2*pi*y);
end

% plot the initial conditions
subplot(1,2,1)
[X,Y]=get_XY(cells);
plot_uns(X,Y,u);

ec=0;
edges=[];
for i=1:3:3*length(tria)
    tidx=(i-1)/3+1;
    for k=0:2
        ss.x=vert(tria(tidx,mod(k  ,3)+1),1);
        ss.y=vert(tria(tidx,mod(k  ,3)+1),2);
        se.x=vert(tria(tidx,mod(k+1,3)+1),1);
        se.y=vert(tria(tidx,mod(k+1,3)+1),2);
        if( new_edge(edges,ec,ss,se) )
            ec=ec+1;
            edges(ec).id=ec;
            edges(ec).ss=ss;
            edges(ec).se=se;
            edges(ec).cen.x=0.5*(edges(ec).ss.x+edges(ec).se.x);
            edges(ec).cen.y=0.5*(edges(ec).ss.y+edges(ec).se.y);
            [nid,pid]=find_nbr(edges(ec).ss,edges(ec).se,ec);
            edges(ec).n=nid;
            edges(ec).p=pid;
        end
    end
end

% for e = edges
%     fprintf('id %d, cen (%.2f,%.2f), pid %d, nid %d \n',e.id,e.cen.x,e.cen.y,e.p,e.n);
% end
% 
% for c = cells
%     fprintf('id %d, f1 %d, f2 %d, f3 %d \n',c.id,c.faces(1).id,c.faces(2).id,c.faces(3).id);
% end

ue=zeros(length(edges),1);
for e = edges
    % get the cells
    if(e.n==0 || e.p==0)
        if(e.n>0)
            ue(e.id)=u(e.n);
        else
            ue(e.id)=u(e.p);
        end
    else
        p = cells(e.p);
        n = cells(e.n);

        np = get_length(n.cen,p.cen);
        nf = get_length(n.cen,e.cen);    
        alpha=abs(nf)/abs(np);
        ue(e.id)=alpha*u(e.p)+(1-alpha)*u(e.n);
    end
    
end

% calculate the gradient at all 

% plot the initial conditions
subplot(1,2,2)
[X,Y]=get_XY(edges);
plot_uns(X,Y,ue);


%% FUNCTIONS

% calculates the distance between two points
% used for side length
function len = get_length(v1,v2)
    len = sqrt( (v1.x-v2.x)^2 + (v1.y-v2.y)^2 );
end

% calculates the cell volume from vertices
function vol = get_volume(v1,v2,v3)
    vol = 0.5*abs( v1.x*(v2.y-v3.y) ... 
                 + v2.x*(v3.y-v1.y) ...
                 + v3.x*(v1.y-v2.y) );
end

% calculates the cell centroid from vertices
function cen = get_centroid(v1,v2,v3)
    cen.x=(1/3)*(v1.x+v2.x+v3.x);
    cen.y=(1/3)*(v1.y+v2.y+v3.y);
end

% get vectors of the centroids
function [X,Y]=get_XY(cells)
 num=length(cells);
 X=zeros(num,1);
 Y=zeros(num,1);
 for c = cells
     X(c.id)=c.cen.x;
     Y(c.id)=c.cen.y;
 end
end

% plot the result
function plot_uns(X,Y,u)
    tri=delaunay(X,Y);
    [r,c]=size(tri);
    h=trisurf(tri,X,Y,u);
%     shading interp;
    axis vis3d;
end

% find the cells that have this edge
function [nid,pid] = find_nbr(es,ee,i)
    global cells
    nid=0;
    pid=0;
    % for every cell    
    for c = cells
        % for every face of every cell
        for f = c.faces
            bool = check_face_eq(f.ss,f.se,es,ee);
            if(bool && nid==0)
                nid=c.id; % this edge has cell id as a nbr
                cells(c.id).faces(f.num).id=i;   % this face corresponds to edge i
            elseif(bool) % if nid has already been assigned...
                pid=c.id;
                cells(c.id).faces(f.num).id=i;
            end
            if(pid ~= 0 && nid ~=0)
                return;
            end
        end
    end
end

% check whether a face and edge are the same
function bool = check_face_eq(fss,fse,es,ee)
    bool=0;
    % check if one of the face vertices is one of the edge vertices
    if(check_vert_eq(fss,es) || check_vert_eq(fss,ee))
        % check if the other face vertex is one of the edge vertices
        if(check_vert_eq(fse,es) || check_vert_eq(fse,ee))
            bool=1; % the edge and the face are the same
        end
    end
end

% check if two points (vertices) are equivalent
function bool = check_vert_eq(v1,v2)
    bool=0;
    if(v1.x==v2.x && v1.y==v2.y)
        bool=1;
    end
end

% check if ss,se is in the list of edges
function bool = new_edge(edges,ec,ss,se)
    bool=1;
    % loop through the known edges
    for i=1:ec
        % check equivalence
        if(check_face_eq(edges(i).ss,edges(i).se,ss,se))
            bool=0;
            return;
        end
    end
end

function plot_mesh

    % plot the unstructured mesh
    figure;
    patch('faces',tria(:,1:3),'vertices',vert, ...
        'facecolor','w', ...
        'edgecolor',[.2,.2,.2]) ;
    hold on; axis image off;
    patch('faces',edge(:,1:2),'vertices',node, ...
        'facecolor','w', ...
        'edgecolor',[.1,.1,.1], ...
        'linewidth',1.5) ;

    for c = cells
        cell_num=sprintf('%d',c.id);
        tc(c.id) = text(c.cen.x,c.cen.y,cell_num);
        tc(c.id).Color='red';
        tc(c.id).FontSize=14;
    end

    for e = edges
        edge_num=sprintf('%d',e.id);
        te(e.id) = text(e.cen.x,e.cen.y,edge_num);
        te(e.id).Color='blue';
    end

    drawnow;

    set(figure(1),'units','normalized', ...
        'position',[.05,.50,.30,.35]) ;

end


