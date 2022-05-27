function plot_mesh_normals( vert,etri,tria,tnum,...
                            node,edge,edges,cells,...
                            hfun)

% plot the mesh with cell ids and visualize face normal fectors
  
% plot the unstructured mesh

figure;
hold on
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

hold on;
for c = cells
    for f = c.faces
        if(strcmp(f.bd,'NA'))
            xbegin=f.cen.x;
            ybegin=f.cen.y;
            plot(xbegin,ybegin,'ko','MarkerSize',3,'MarkerFaceColor','k');
            dx=min(abs(xbegin),min(abs(1-xbegin),hfun/10));
            dy=min(abs(ybegin),min(abs(1-ybegin),hfun/10));
            xend=xbegin+dx*f.nv.x;
            yend=ybegin+dy*f.nv.y;
            plot([xbegin,xend],[ybegin,yend],'k-');
        end
    end
end


end

