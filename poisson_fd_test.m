% Finite Differece code to solve poisson equation
% 7/4/18
% Nick Stegmeier
clear all; close all; clc;

idx=1;
for N=100:10:100

% Nx = 40;
% Ny = 40;
Nx=N;
Ny=N;

dx=1/Nx;
dy=1/Ny;

x=dx/2:dx:1-dx/2;
y=dy/2:dy:1-dy/2;

g=zeros(Nx,Ny);
for j=1:Ny
    for i=1:Nx
%         g = @(x,y) 20*cos(3*pi*x)*sin(2*pi*y);
        g(j,i)=20*cos(pi*(x(i)+1/2))*sin(2*pi*y(j));
%           g(i,j)=20*cos(3*pi*x(i))*sin(2*pi*y(i));
%           g(i,j) = -1.0;
    end
end

% u on xfaces
for j=1:Ny
    uxf(1,j)=0.0;
    uxf(2,j)=0.0;
end

% u on yfaces
for i=1:Nx
    uyf(i,1)=0.0;
    uyf(i,2)=0.0;
end

% figure(1);
[X,Y]=meshgrid(x,y);
% surf(X,Y,g);

A=zeros(Nx*Ny,Nx*Ny);
b=zeros(Nx*Ny,1);
for j=1:Ny
    for i=1:Nx
        
        ap=-2*(1/dx^2 + 1/dy^2);
        aex=1/dx^2;
        awx=1/dx^2;
        aey=1/dy^2;
        awy=1/dy^2;
        
        aexf=1;
        awxf=1;
        aeyf=1;
        awyf=1;
        
        if(i==1)
            awxf=0;
        end
        if(i==Nx)
            aexf=0;
        end
        if(j==1)
            awyf=0;
        end
        if(j==Ny)
            aeyf=0;
        end       
        
        row=i+Nx*(j-1);
        col=i+Nx*(j-1);
        A(row,col)=ap;
        b(row)=g(Ny-j+1,i);
        if(awxf) 
            A(row,col-1)=awx;
        else
            b(row)=b(row)-2*uxf(1,j);
        end
        if(aexf) 
            A(row,col+1)=aex;
        else
            b(row)=b(row)-2*uxf(2,j);
        end
        if(awyf) 
            A(row,col-Nx)=awy;
        else
            b(row)=b(row)-2*uyf(i,1);
        end
        if(aeyf) 
            A(row,col+Ny)=aey;
        else
            b(row)=b(row)-2*uyf(i,2);
        end
        
    end
end

u=zeros(Nx*Ny,1);
u = A\b;

for j=1:Ny
    for i=1:Nx
        row=i+Nx*(j-1);
        ud(Ny-(j-1),Nx-(i-1))=u(row);
    end
end

figure(2);
surf(X,Y,ud)
N
umaxs(idx) = max(max(u));
idx=idx+1;

end