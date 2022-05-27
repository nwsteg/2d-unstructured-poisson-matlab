function [edges,cells] = initialize_edges( vert,etri,tria,tnum,cells )

ec=0;
edges=[];
for i=1:3:3*length(tria)
    tidx=(i-1)/3+1;
    for k=0:2
        ss.x=vert(tria(tidx,mod(k  ,3)+1),1);
        ss.y=vert(tria(tidx,mod(k  ,3)+1),2);
        se.x=vert(tria(tidx,mod(k+1,3)+1),1);
        se.y=vert(tria(tidx,mod(k+1,3)+1),2);
        if( new_edge(edges,ec,ss,se) )
            ec=ec+1;
            fprintf('\tec=%d\n',ec);
            edges(ec).id=ec;
            edges(ec).ss=ss;
            edges(ec).se=se;
            edges(ec).cen.x=0.5*(edges(ec).ss.x+edges(ec).se.x);
            edges(ec).cen.y=0.5*(edges(ec).ss.y+edges(ec).se.y);
            [nid,pid,cells]=find_nbr(edges(ec).ss,edges(ec).se,ec,cells);
            edges(ec).n=nid;
            edges(ec).p=pid;
        end
    end
end

end

