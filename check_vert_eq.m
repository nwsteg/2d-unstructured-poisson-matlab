% check if two points (vertices) are equivalent
function bool = check_vert_eq(v1,v2)
    bool=0;
    if(v1.x==v2.x && v1.y==v2.y)
        bool=1;
    end
end