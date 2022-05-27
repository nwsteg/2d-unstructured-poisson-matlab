function [bdry,bool] = check_boundary_face( face )

bdry = 'NA';
bool = 0;

% check if a particle face is entirely on the boundary
% if yes, return which boundary (east, west, north, south)
% else, return NA

% west boundary
if( face.ss.x == 0 && face.se.x == 0)
    bdry = 'west';
    bool = 1;
    return;
end

% bottom boundary
if( face.ss.y == 0 && face.se.y == 0)
    bdry = 'south';
    bool = 1;
    return;
end

% right boundary
if( face.ss.x == 1 && face.se.x == 1)
    bdry = 'east';
    bool = 1;
    return;
end

% top boundary
if( face.ss.y == 1 && face.se.y == 1)
    bdry = 'north';
    bool = 1;
    return;
end;

return;

end

