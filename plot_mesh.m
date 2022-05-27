function plot_mesh(vert,etri,tria,tnum,node,edge,edges,cells)

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
% 
    for e = edges
        edge_num=sprintf('%d',e.id);
        te(e.id) = text(e.cen.x,e.cen.y,edge_num);
        te(e.id).Color='blue';
    end

    drawnow;

    set(figure(1),'units','normalized', ...
        'position',[.05,.50,.30,.35]) ;

end