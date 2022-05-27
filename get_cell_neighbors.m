function cells = get_cell_neighbors(cells,edges)
% this function finds the neighbors of each cell
% it also determines if a given cell is "near" the domain boundary
% ie, if a cell touches another cell that is on the boundary

% determine the cell neighor normal to each face
% loop through all edges
for e = edges
    flag=0;
    % get cells touching edge e
    % if either e.p or e.n is zero, the edge is on the boundary
    % thus it doesn't make sense to identify the touching cells as
    % neighbors
    if(e.p~=0); p = cells(e.p); else; continue; end;
    if(e.n~=0); n = cells(e.n); else; continue; end;
    %fprintf('edge %d\n',e.id)
    % loop over faces of p
    for fp = p.faces
        % loop over faces of n
        for fn = n.faces
%             fprintf('\t fp (%d,%d), fn (%d,%d) \n',fp.id,fp.num,fn.id,fn.num)
            % check if faces are the same
            if(fp.id==fn.id)
%                 fprintf('%d and %d are neighbors \n',n.id,p.id)
                % assign e.p.face to have nbr p.id
                cells(e.p).faces(fp.num).nbr=n.id;
                % assign e.n.face to have nbr n.id
                cells(e.n).faces(fn.num).nbr=p.id;
                % turn on flag so that we exit to edges loop
                flag=1;
                break;
            end
        end
        if(flag==1)
            break;
        end
    end
end

% now get a list of all neighbors for each cell
for c = cells
    
    c.is_near_bd=0;
    c.bd=0;

    for f = c.faces
        if( ~strcmp(f.bd,'NA') )
            c.bd=1;
        end
        c.nbrs=[c.nbrs f.nbr];
    end
    
    % check the neighbors and see if they're on the boundary
    for c_nbr = cells(find(c.nbrs>0))
        if(c_nbr.isbd==1)
            c.is_near_bd=1;
        end
    end
     
    cells(c.id).nbrs=c.nbrs;
    cells(c.id).bd=c.bd;
    cells(c.id).is_near_bd=c.is_near_bd;
    
    % assign the is_near_bd flag for each face of cell c
    for f = c.faces
        flag=0;
        if(f.nbr ~= 0)
            if( cells(f.nbr).isbd==1 )
                flag=1;
            end
        end
        cells(c.id).faces(f.num).is_near_bd=flag;
    end
    
end

end

