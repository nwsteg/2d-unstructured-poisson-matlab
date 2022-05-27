function print_cells( cells )

for c = cells
    fprintf('\ncell id = %3d at (%.2f,%.2f) with isbd=%d\n',c.id,c.cen.x,c.cen.y,c.isbd);
    fprintf('\t neighbors = ');
    for nbr = c.nbrs
        fprintf('%d ',nbr)
    end
    fprintf('\n')
    for f = c.faces
        fprintf('\t face id = %3d, center (%.2f,%.2f), boundary %8s,     near_boundary %d, nbr id %4d, alpha %f \n',f.id,f.cen.x,f.cen.y,f.bd,f.is_near_bd,f.nbr,f.alpha);
    end
end

end

