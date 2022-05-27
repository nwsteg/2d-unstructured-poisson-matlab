function bool = new_edge(edges,ec,ss,se)
% check if ss,se is in the list of edges
    bool=1;
    % loop through the known edges
    for i=1:ec
        % check equivalence
        if(check_face_eq(edges(i).ss,edges(i).se,ss,se))
            bool=0;
            return;
        end
    end
end
