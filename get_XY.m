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