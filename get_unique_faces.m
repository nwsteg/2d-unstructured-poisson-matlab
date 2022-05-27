% Finds the face indices s.t. u1.faces(fnum1) and u1.faces(fnum2)
% are distinct from face u0_id

% INPUTS
%   u0_id the face id from cell u0 to check against
%   u1 the cell adjacent to u0 which we are comparing

% OUTPUS
%   fnum1 the index in u1.faces s.t. u1.faces(fnum1).id ~= u0_id
%   fnum2 the index in u2.faces s.t. fnum2 ~= fnum1 and u1.faces(fnum2).id
%   ~= u0_id

function [fnum1,fnum2] = get_unique_faces(u0_id,u1)
    u1_ids(1)=u1.faces(1).id;
    u1_ids(2)=u1.faces(2).id;
    u1_ids(3)=u1.faces(3).id;
    fnums_out(1:2)=0;
    idx=1;
    num=1;
    for fid = u1_ids
        if(fid ~= u0_id)
            fnums_out(idx)=num;
            idx=idx+1;
        end
        num=num+1;
    end
    fnum1=fnums_out(1);
    fnum2=fnums_out(2);
end