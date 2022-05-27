function cells = calculate_alphas( cells )
% calculate_alphas -- find the proportion-distance for each cell face
%   To calculate the value at a face, you take the "average"
%   value between the cells that share the face
%   This average value, alpha, is calculated by finding the ratio of
%   the distance from each cell.
%   The formula is explained on CFD wiki and in the pdf.

for P = cells
    for f = P.faces
        if( f.nbr ~=0 )
            N = cells(f.nbr);
            nbr_to_face = get_length(N.cen,f.cen);
            total = get_length(P.cen,N.cen);

            % assign alpha for each cell face
            cells(P.id).faces(f.num).alpha = abs(nbr_to_face/total);
        end
    end
end


end

