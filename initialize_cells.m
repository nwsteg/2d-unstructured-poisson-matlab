function cells = initialize_cells( vert,etri,tria,tnum )

num_cells=length(tria);

% initialize mesh structure
cells = []; % cell
cells.vert=[]; % cell vertices
cells.faces=[];    % cell faces
cells.faces.ss=[]; % cell face start
cells.faces.se=[]; % cell face end
cells.faces.l=0.0; % side length
cells.vol=0.0; % cell volume
cells.cen=[]; % cell centroid
cells.isbd=0;
cells.nbrs=[];


% populate the mesh data structure
for i=1:num_cells

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
    
    % general info
    cells(i).faces(1).num=1;
    cells(i).faces(1).id=0;
    cells(i).faces(1).nbr=0;
    cells(i).faces(1).ss=cells(i).vert(1);
    cells(i).faces(1).ss.id=cells(i).vert(1).id;
    cells(i).faces(1).se=cells(i).vert(2);
    cells(i).faces(1).se.id=cells(i).vert(2).id;
    cells(i).faces(1).cen=get_midpoint(cells(i).faces(1).ss,cells(i).faces(1).se);
    cells(i).faces(1).l=get_length(cells(i).faces(1).ss,...
                                   cells(i).faces(1).se); 
                               
    % face normal
    delx=cells(i).faces(1).ss.x-cells(i).faces(1).se.x;
    dely=cells(i).faces(1).ss.y-cells(i).faces(1).se.y;
    alpha=sqrt(delx*delx+dely*dely);
    cells(i).faces(1).nv.x =  dely/alpha;
    cells(i).faces(1).nv.y = -delx/alpha;
    
    % boundary?
    [cells(i).faces(1).bd, cells(i).faces(1).isbd] = check_boundary_face( cells(i).faces(1) );
                   
    % general info
    cells(i).faces(2).num=2;
    cells(i).faces(2).id=0;
    cells(i).faces(2).nbr=0;
    cells(i).faces(2).ss=cells(i).vert(2);
    cells(i).faces(2).ss.id=cells(i).vert(2).id;
    cells(i).faces(2).se=cells(i).vert(3);
    cells(i).faces(2).se.id=cells(i).vert(3).id;
    cells(i).faces(2).cen=get_midpoint(cells(i).faces(2).ss,cells(i).faces(2).se);
    cells(i).faces(2).l=get_length(cells(i).faces(2).ss,...
                                   cells(i).faces(2).se);
    % face normal
    delx=cells(i).faces(2).ss.x-cells(i).faces(2).se.x;
    dely=cells(i).faces(2).ss.y-cells(i).faces(2).se.y;
    alpha=sqrt(delx*delx+dely*dely);
    cells(i).faces(2).nv.x =  dely/alpha;
    cells(i).faces(2).nv.y = -delx/alpha;
    
    % boundary?
    [cells(i).faces(2).bd, cells(i).faces(2).isbd] = check_boundary_face( cells(i).faces(2) );
    
    % general info
    cells(i).faces(3).num=3; 
    cells(i).faces(3).id=0;
    cells(i).faces(3).nbr=0;
    cells(i).faces(3).ss=cells(i).vert(3);
    cells(i).faces(3).ss.id=cells(i).vert(3).id;
    cells(i).faces(3).se=cells(i).vert(1);
    cells(i).faces(3).se.id=cells(i).vert(1).id;
    cells(i).faces(3).cen=get_midpoint(cells(i).faces(3).ss,cells(i).faces(3).se);
    cells(i).faces(3).l=get_length(cells(i).faces(3).ss,...
                                   cells(i).faces(3).se); 
                               
    % face normal
    delx=cells(i).faces(3).ss.x-cells(i).faces(3).se.x;
    dely=cells(i).faces(3).ss.y-cells(i).faces(3).se.y;
    alpha=sqrt(delx*delx+dely*dely);
    cells(i).faces(3).nv.x =  dely/alpha;
    cells(i).faces(3).nv.y = -delx/alpha;
    
    % boundary?
    [cells(i).faces(3).bd, cells(i).faces(3).isbd] = check_boundary_face( cells(i).faces(3) );
    
                           
    % calculate the volume
    cells(i).vol=get_volume(cells(i).vert(1), ...
                            cells(i).vert(2), ...
                            cells(i).vert(3) );
                        
    % calculate the centroid
    cells(i).cen=get_centroid(cells(i).vert(1), ...
                              cells(i).vert(2), ...
                              cells(i).vert(3) );
           
    cells(i).isbd=0;
    for fc=cells(i).faces
        if ( ~strcmp(fc.bd,'NA') )
            cells(i).isbd=1;
        end
    end
    
end


end

