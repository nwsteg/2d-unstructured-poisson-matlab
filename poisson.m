% 2D Poisson Equation
% Unstructured Mesh
% Finite Volume Discretization

% Nick Stegmeier
% 7/3/2018

clear all; close all; clc;

zero.x=0;
zero.y=0;
tic
% make the mesh
hfun = 0.05; % "target" edge length (1x1 domain)
plot = 1; % plot the wireframe (yes=1)
[vert,etri,tria,tnum,node,edge] = generate_mesh(hfun,plot);
fprintf('There are %d cells and %d edges... \n',length(tria),length(etri))

% calculate mesh quantities (id, centroid, faces, etc)
fprintf('Initializing cells...\n')
cells = initialize_cells(vert,etri,tria,tnum );

% identify and label all edges
% label faces in cells data structure
fprintf('Initializing edges... \n')
[edges,cells] = initialize_edges( vert,etri,tria,tnum,cells );
%print_edges(edges);

% determine the cell neighor normal to each face
fprintf('Finding neighbors... \n')
cells = get_cell_neighbors(cells,edges);

% get geometric weighting factor for each cell face
fprintf('Calculating alphas... \n')
cells = calculate_alphas(cells);

% check basic cell data
print_cells(cells);

% look at the mesh with cell and face ids
% figure(1)
% plot_mesh(vert,etri,tria,tnum,node,edge,edges,cells);

% check the cell face normals 
% plot_mesh_normals(vert,etri,tria,tnum,node,edge,edges,cells,hfun);
num_cells=length(cells);

toc
if 1==1
g = @(x,y) 20*cos(pi*(x+1/2))*sin(2*pi*y);
% g = @(x,y) 20*cos(3*pi*x)*sin(2*pi*y);
gv = zeros(num_cells,1);
ubd=0; % boundary condition for u
for c = cells
    x=c.cen.x;
    y=c.cen.y;
    gv(c.id) = g(x,y);
end

[X,Y]=get_XY(cells);
%figure(2)
%plot_uns(X,Y,gv);

fprintf('Constructing matrix... \n')
A=zeros(num_cells,num_cells);
b=zeros(num_cells,1);
total_vol=0;
for c = cells
    
    % reset cells
    u0=[]; u1=[]; u2=[]; u3=[];
    u4=[]; u5=[]; u6=[]; u7=[];
    u8=[]; u9=[];
    
    % reset faces
    f1=[]; f2=[]; f3=[]; f4=[];
    f5=[]; f6=[]; f7=[]; f8=[];
    f9=[];
    
    % reset normals
    P1=[]; P2=[]; P3=[];
    N4=[]; N5=[]; N6=[];
    N7=[]; N8=[]; N9=[];
    
    % reset rhs face vals
    
    % matrix coefficient flags
    Cu0f=1;
    Cu1f=0; Cu4f=0; Cu5f=0;
    Cu2f=0; Cu6f=0; Cu7f=0;
    Cu3f=0; Cu8f=0; Cu9f=0;
    
    % matrix & rhs coefficients
    Cu0=0;
    Cu1=0; Cu4=0; Cu5=0;
    Cu2=0; Cu6=0; Cu7=0;
    Cu3=0; Cu8=0; Cu9=0;
    rhs=0;
    
    % beta values -- geometric interaction coefficients
    b12=0; b13=0; b14=0; b15=0;
    b21=0; b23=0; b26=0; b27=0;
    b31=0; b32=0; b38=0; b39=0;
    
    % alpha values -- geometric proportions
    alpha1=0; alpha4=0; alpha5=0;
    alpha2=0; alpha6=0; alpha7=0;
    alpha3=0; alpha8=0; alpha9=0;
    
    % main cell
    u0=c;
    rhs = rhs + u0.vol*g(u0.cen.x,u0.cen.y);
    
    % faces
    f1 = u0.faces(1);
    f2 = u0.faces(2);
    f3 = u0.faces(3);
    
    % face area normals = area*vector
    P1 = sm(f1.l,f1.nv);
    P2 = sm(f2.l,f2.nv);
    P3 = sm(f3.l,f3.nv);
    
    % ----- face 1 -----
    if f1.isbd==true
        vec.x = u0.cen.x - f1.cen.x;
        vec.y = u0.cen.y - f1.cen.y;
        dz = get_length(u0.cen,f1.cen);     
%         uf1 = g(f1.cen.x,f1.cen.y);
        uf1=ubd;
        Cu0 = Cu0 + vec_dot(vec,P1)/(dz*dz);
        rhs = rhs + uf1*vec_dot(vec,P1)/(dz*dz);
    else
        Cu1f=1;
        u1 = cells(u0.nbrs(1));
        omega1 = u0.vol + u1.vol;
        alpha1 = f1.alpha;
        [fnum1,fnum2] = get_unique_faces(f1.id,u1);
        f4 = u1.faces(fnum1); 
        f5 = u1.faces(fnum2); 
        N4 = sm(f4.l,f4.nv);
        N5 = sm(f5.l,f5.nv);
        b12 = (1/omega1)*vec_dot(P1,P2);
        b13 = (1/omega1)*vec_dot(P1,P3);
        b14 = (1/omega1)*vec_dot(P1,N4);
        b15 = (1/omega1)*vec_dot(P1,N5);
        if f4.isbd == true
%             uf4 = g(f4.cen.x,f4.cen.y);
            uf4 = ubd;
            rhs = rhs - b14*uf4;
            b14 = 0;
        else
            Cu4f=1;
            alpha4 = f4.alpha;
            u4 = cells(f4.nbr);
        end
        if f5.isbd == true
%             uf5 = g(f5.cen.x,f5.cen.y);
            uf5 = ubd;
            rhs = rhs - b15*uf5;
            b15 = 0;
        else
            Cu5f=1;
            alpha5 = f5.alpha;
            u5 = cells(f5.nbr);
        end
    end
    % -----
    
    % ----- face 2 -----
    if f2.isbd==true
        vec.x = u0.cen.x - f2.cen.x;
        vec.y = u0.cen.y - f2.cen.y;
        dz = get_length(u0.cen,f2.cen);     
%         uf2 = g(f2.cen.x,f2.cen.y);
        uf2 = ubd;
        Cu0 = Cu0 + vec_dot(vec,P2)/(dz*dz);
        rhs = rhs + uf2*vec_dot(vec,P2)/(dz*dz);
    else
        Cu2f=1;
        u2 = cells(u0.nbrs(2));
        omega2 = u0.vol + u2.vol;
        alpha2 = f2.alpha;
        [fnum1,fnum2] = get_unique_faces(f2.id,u2);
        f6 = u2.faces(fnum1);
        f7 = u2.faces(fnum2);
        N6 = sm(f6.l,f6.nv);
        N7 = sm(f7.l,f7.nv);
        b21 = (1/omega2)*vec_dot(P2,P1);
        b23 = (1/omega2)*vec_dot(P2,P3);
        b26 = (1/omega2)*vec_dot(P2,N6);
        b27 = (1/omega2)*vec_dot(P2,N7);
        if f6.isbd == true
%             uf6 = g(f6.cen.x,f6.cen.y);
            uf6 = ubd;
            rhs = rhs - b26*uf6;
            b26 = 0;
        else
            Cu6f=1;
            alpha6 = f6.alpha;
            u6 = cells(f6.nbr);
        end
        if f7.isbd == true
%             uf7 = g(f7.cen.x,f7.cen.y);
            uf7 = ubd;
            rhs = rhs - b27*uf7;
            b27 = 0;
        else
            Cu7f=1;
            alpha7 = f7.alpha;
            u7 = cells(f7.nbr);
        end    
    end
    
    % ----- face 3 -----
    if f3.isbd == true
        vec.x = u0.cen.x - f3.cen.x;
        vec.y = u0.cen.y - f3.cen.y;
        dz = get_length(u0.cen,f3.cen);     
%         uf3 = g(f3.cen.x,f3.cen.y);
        uf3 = ubd;
        Cu0 = Cu0 + vec_dot(vec,P3)/(dz*dz);
        rhs = rhs + uf3*vec_dot(vec,P3)/(dz*dz);
    else
        Cu3f=1;
        u3 = cells(u0.nbrs(3));
        omega3 = u0.vol + u3.vol;
        alpha3 = f3.alpha;
        [fnum1,fnum2] = get_unique_faces(f3.id,u3);
        f8 = u3.faces(fnum1);
        f9 = u3.faces(fnum2);
        N8 = sm(f8.l,f8.nv);
        N9 = sm(f9.l,f9.nv);
        b31 = (1/omega3)*vec_dot(P3,P1);
        b32 = (1/omega3)*vec_dot(P3,P2);
        b38 = (1/omega3)*vec_dot(P3,N8);
        b39 = (1/omega3)*vec_dot(P3,N9);
        if f8.isbd == true
%             uf8 = g(f8.cen.x,f8.cen.y);
            uf8 = ubd;
            rhs = rhs - b38*uf8;
            b38 = 0;
        else
            Cu8f=1;
            alpha8 = f8.alpha;
            u8 = cells(f8.nbr);
        end
        if f9.isbd == true
%             uf9 = g(f9.cen.x,f9.cen.y);
            uf9 = ubd;
            rhs = rhs - b39*uf9;
            b39 = 0;
        else
            Cu9f=1;
            alpha9 = f9.alpha;
            u9 = cells(f9.nbr);
        end
    end
    
    Cu0=Cu0+(b21+b31)*alpha1+(b12+b32)*alpha2+(b13+b23)*alpha3;
    Cu1=Cu1+(b21+b31)*(1-alpha1)+b14*alpha4+b15*alpha5;
    Cu2=Cu2+(b12+b32)*(1-alpha2)+b26*alpha6+b27*alpha7;
    Cu3=Cu3+(b13+b23)*(1-alpha3)+b38*alpha8+b39*alpha9;
    Cu4=b14*(1-alpha4);
    Cu5=b15*(1-alpha5);
    Cu6=b26*(1-alpha6);
    Cu7=b27*(1-alpha7);
    Cu8=b38*(1-alpha8);
    Cu9=b39*(1-alpha9);
    
    A(u0.id,u0.id)=Cu0;
    if(Cu1f == true); A(u0.id,u1.id)=Cu1; end;
    if(Cu2f == true); A(u0.id,u2.id)=Cu2; end;
    if(Cu3f == true); A(u0.id,u3.id)=Cu3; end;
    if(Cu4f == true); A(u0.id,u4.id)=Cu4; end;
    if(Cu5f == true); A(u0.id,u5.id)=Cu5; end;
    if(Cu6f == true); A(u0.id,u6.id)=Cu6; end;
    if(Cu7f == true); A(u0.id,u7.id)=Cu7; end;
    if(Cu8f == true); A(u0.id,u8.id)=Cu8; end;
    if(Cu9f == true); A(u0.id,u9.id)=Cu9; end;
    
    b(u0.id) = rhs;
    
end

sol = A\b;
figure(2)
plot_uns(X,Y,sol);
end
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%