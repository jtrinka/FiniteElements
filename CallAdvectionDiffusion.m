%% CallAdvectionDiffusion.m code
clear; clc;
%This code loads a domain that is created by the user from using the
%ChooseDomain.m and DomainConstruction.m files. It then runs the domain
%through the AdvectionDiffusionFiniteElements.m code. Use the load command to
%load your data coordinates data, nodes data, Neumann data, Dirichlet data, minimum bound on the domain, and
%maximum bound on the domain.
load xyfirstmesh.mat
load nodesfirstmesh.mat
load neumannfirstmesh.mat
load Dirichfirstmesh.mat
load savexsizeminfirstmesh.mat
load saveysizeminfirstmesh.mat
load savexsizemaxfirstmesh.mat
load saveysizemaxfirstmesh.mat
%Initial Condition, velocity vector field, diffusivity term, maximum time
%to run, and time step
constant = (1/(sqrt(2*pi)^2*1));
f=@(x,y) constant*exp(-.0007*(x-75).^2-.0007*(y-155).^2);
v=[30;30];
D=1000;
T=20; dt=0.01;
AdvectionDiffusionWithNeumannWorking(xy,nodes,neumann,Dirich,savexsizemin,saveysizemin,savexsizemax,saveysizemax,f,v,D,T,dt)