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

