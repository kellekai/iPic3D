clear all;
close all

Lx=20;nxc=128;dx=Lx/(nxc+1);
Ly=20;nyc=128;dy=Ly/(nyc+1);
Lz=1;nzc=1;dz=Lz/(nzc+1);

x = linspace(0,Lx,nxc+1);
y = linspace(0,Ly,nyc+1);
z = linspace(0,Lz,nzc+1);

[X,Y,Z] = ndgrid(x,y,z);
[Nx,Ny,Nz]=size(X)
% 
% X=permute(X,[3 2 1]);
% Y=permute(Y,[3 2 1]);
% Z=permute(Z,[3 2 1]);


ns=2;
B0=1; L=0.1; 

N0=[1;1]/4/pi;
qom=[-256;1];

U0=[0;0];V0=[0;0];W0=[.00325;-.01624];


Bx=B0*tanh((Y-Ly/2)./L);
By=zeros(Nx,Ny,Nz);
Bz=zeros(Nx,Ny,Nz);

Ex=zeros(Nx,Ny,Nz);
Ey=zeros(Nx,Ny,Nz);
Ez=zeros(Nx,Ny,Nz);

rho = sech((Y-Ly/2)./L).^2;


!rm Initial-Fields_000000.h5
opath='Initial-Fields_000000.h5'
h5create(opath,'/Step#0/Block/Bx/0',[Nx, Ny, Nz]);
h5write(opath,'/Step#0/Block/Bx/0',Bx)

h5create(opath,'/Step#0/Block/By/0',[Nx, Ny, Nz]);
h5write(opath,'/Step#0/Block/By/0',By)

h5create(opath,'/Step#0/Block/Bz/0',[Nx, Ny, Nz]);
h5write(opath,'/Step#0/Block/Bz/0',Bz)

h5create(opath,'/Step#0/Block/Ex/0',[Nx, Ny, Nz]);
h5write(opath,'/Step#0/Block/Ex/0',Ex)

h5create(opath,'/Step#0/Block/Ey/0',[Nx, Ny, Nz]);
h5write(opath,'/Step#0/Block/Ey/0',Ey)

h5create(opath,'/Step#0/Block/Ez/0',[Nx, Ny, Nz]);
h5write(opath,'/Step#0/Block/Ez/0',Ez)

for is=1:ns
h5create(opath,['/Step#0/Block/Jx_' num2str(is-1) '/0'],[Nx, Ny, Nz]);
h5write(opath,['/Step#0/Block/Jx_' num2str(is-1) '/0'],N0(is)*rho*U0(is)*qom(is))

h5create(opath,['/Step#0/Block/Jy_' num2str(is-1) '/0'],[Nx, Ny, Nz]);
h5write(opath,['/Step#0/Block/Jy_' num2str(is-1) '/0'],N0(is)*rho*V0(is)*qom(is))

h5create(opath,['/Step#0/Block/Jz_' num2str(is-1) '/0'],[Nx, Ny, Nz]);
h5write(opath,['/Step#0/Block/Jz_' num2str(is-1) '/0'],N0(is)*rho*W0(is)*qom(is))

h5create(opath,['/Step#0/Block/rho_' num2str(is-1) '/0'],[Nx, Ny, Nz]);
h5write(opath,['/Step#0/Block/rho_' num2str(is-1) '/0'],N0(is)*rho*qom(is))
end
h5writeatt(opath,'/Step#0','nspec',int32(ns));