% calculates the cell centroid from vertices
function cen = get_centroid(v1,v2,v3)
    cen.x=(1/3)*(v1.x+v2.x+v3.x);
    cen.y=(1/3)*(v1.y+v2.y+v3.y);
end
