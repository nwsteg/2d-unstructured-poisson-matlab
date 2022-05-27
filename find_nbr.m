function [nid,pid,cells] = find_nbr(es,ee,i,cells)
    nid=0;
    pid=0;
    % for every cell    
    for c = cells
        % for every face of every cell
        for f = c.faces
            bool = check_face_eq(f.ss,f.se,es,ee);
            if(bool && nid==0)
                nid=c.id; % this edge has cell id as a nbr
                cells(c.id).faces(f.num).id=i;   % this face corresponds to edge i
            elseif(bool) % if nid has already been assigned...
                pid=c.id;
                cells(c.id).faces(f.num).id=i;
            end
            if(pid ~= 0 && nid ~=0)
                return;
            end
        end
    end
end