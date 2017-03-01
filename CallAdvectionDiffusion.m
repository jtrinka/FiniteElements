%% CallAdvectionDiffusion.m code
clear; clc;
%This code loads a domain that is created by the user from using the
%ChooseDomain.m and DomainConstruction.m files. It then runs the domain
%through the AdvectionDiffusionFiniteElements.m code. In order to get access to
%the AdvectionDiffusionFiniteElements.m file, please email me at 
%jordantrinka4@hotmail.com. The AdvectionDiffusionFiniteElements.m file
%contains diffusion code drawn directly from
%https://www.math.hu-berlin.de/~cc/cc_homepage/download/1999-AJ_CC_FS-50_Lines_of_Matlab.pdf
%as well as the original research I produced for the advective term.
%The code I built for the advective term was inspired from the
%previously mentioned source. I also had help from Dr. Eric Sullivan of Carroll College in Helena, Montana. Use the load command to
%load your data coordinates data, nodes data, Neumann data, Dirichlet data, minimum bound on the domain, and
%maximum bound on the domain.

%First domain
% load xyfirstmesh.mat
% load nodesfirstmesh.mat
% load neumannfirstmesh.mat
% %load Dirichfirstmesh.mat
% load savexsizeminfirstmesh.mat
% load saveysizeminfirstmesh.mat
% load savexsizemaxfirstmesh.mat
% load saveysizemaxfirstmesh.mat
% load plotufirstmesh.mat
% load plotvfirstmesh.mat
% load myxyfirstmesh.mat

%Second domain
load xysecondmesh.mat
load nodessecondmesh.mat
load neumannsecondmesh.mat
load Dirichsecondmesh.mat
load savexsizeminsecondmesh.mat
load saveysizeminsecondmesh.mat
load savexsizemaxsecondmesh.mat
load saveysizemaxsecondmesh.mat
load plotusecondmesh.mat
load plotvsecondmesh.mat
load myxysecondmesh.mat

neumann1=neumann(1:125,:);
neumann2=neumann(130:end,:);
Dirich1=neumann(126:129,:);
Dirich=[Dirich;Dirich1];
neumann=[neumann1;neumann2];


myvelocities=[reshape(plotu',length(plotu(:,1))^2,1),reshape(plotv',length(plotv(:,1))^2,1)];
nNodes=unique(nodes);
v=zeros(length(nNodes(:,1)),2);
epsilon=5;
for i=1:length(nNodes)
    for j=1:length(myxy(:,1))
        if abs(xy(i,1)-myxy(j,1))<epsilon && abs(xy(i,2)-myxy(j,2))<epsilon
            v(i,1)=myvelocities(j,1);
            v(i,2)=myvelocities(j,2);
        end
        
    end
end
v=v*(1/300);

%Initial Condition, velocity vector field, diffusivity term, maximum time
%to run, and time step
constant = (1/(sqrt(2*pi)^2*170));
f=@(x,y) 0*constant*exp(-.0001*(x-375).^2-.0001*(y-290).^2);
%v=[-20;0];
%try high diffusivity with the neumann boundary condition pulling shit out
%on the left
D=1000;
T=200; dt=.01;
AdvectionDiffusionFiniteElements(xy,nodes,neumann,Dirich,savexsizemin,saveysizemin,savexsizemax,saveysizemax,f,v,D,T,dt)





