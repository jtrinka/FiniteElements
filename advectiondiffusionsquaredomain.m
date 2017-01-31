%% Advection-Diffusion on Square Domain
%All diffusion is drawn directly from
%https://www.math.hu-berlin.de/~cc/cc_homepage/download/1999-AJ_CC_FS-50_Lines_of_Matlab.pdf
%advection code is also inspired by this source.
clear; clc;
x=linspace(0,10,40);
y=linspace(0,10,40);
savexsize = length(x);
saveysize = length(y);
[x,y]=meshgrid(x,y);
x = x(1:end);
y = y(1:end);
x=x';
y=y';
D=.5;
for i=1:1:length(x)
    nodenumber(i,1) = i;
    
end

Data = [nodenumber,x,y];
    
% Triangulate it and assign nodes a number
nodes = delaunay(x,y); %elements3 is nodes
for i=1:1:length(nodes(:,1))
    trinumber(i,1) = i;
    
end
triandnodes=[trinumber,nodes];    
    
%triplot(nodes,x,y) %plot to see it
%text(x,y,num2str(nodenumber(:,1)))
%hold on
for i =1:1:length(nodes(:,1))
   x1=x(nodes(i,1));
   x2=x(nodes(i,2));
   x3=x(nodes(i,3));
   
   y1=y(nodes(i,1));
   y2=y(nodes(i,2));
   y3=y(nodes(i,3));
   
   Bx = (x1+x2+x3)/3;
   By = (y1+y2+y3)/3;
   
   %text(Bx,By,num2str(i))    
end
T=2; dt=0.01; N=T/dt;
v=[0;1];
f=@(x,y) exp(-(x-5.5).^2-(y-2.5).^2);
xy=[x,y]; %coordinates data
%elements3 data is nodes
Dirich = DirichletSteadyState(nodes,0,1,0,1,xy);%Dirichlet Nodes
Freenodes = unique(setdiff(nodes,Dirich)); %Freenodes
A = sparse(size(xy,1),size(xy,1));
B = sparse(size(xy,1),size(xy,1));
C = sparse(size(xy,1),size(xy,1));
U = zeros(size(xy,1),N+1);
%Assembly
for j=1:size(nodes,1)
    A(nodes(j,:),nodes(j,:)) = A(nodes(j,:),nodes(j,:))+stima3(xy(nodes(j,:),:));
end

for j=1:size(nodes,1)
   B(nodes(j,:),nodes(j,:))=B(nodes(j,:),nodes(j,:))+det([1,1,1;xy(nodes(j,:),:)'])*[2,1,1;1,2,1;1,1,2]/24;
end

for j=1:size(nodes,1)
    d = size(xy(nodes(j,:),:),2);
         %gradient of basis functions centered at each of the three nodes
         %in each triangular element
        G = [ones(1,d+1);xy(nodes(j,:),:)']\[zeros(1,d);eye(d)]; %G is fine
    C(nodes(j,:),nodes(j,:))=C(nodes(j,:),nodes(j,:))+det([1,1,1;xy(nodes(j,:),:)'])*reshape([G*v;G*v;G*v],3,3)*(1/6);
end

%Initial Condition
U(:,1) = f(x,y);%zeros(size(xy,1),1);
MyReshapeU = reshape(U(:,1),savexsize,saveysize);
MyReshapeU(1,:)=0;
MyReshapeU(end,:)=0;
MyReshapeU(:,1)=0;
MyReshapeU(:,end)=0;
U(:,1) = reshape(MyReshapeU,savexsize*saveysize,1);
%time steps
 for n=2:N+1
     b=sparse(size(xy,1),1);
%     %volume forces
     for j=1:size(nodes,1)
         b(nodes(j,:)) = b(nodes(j,:));%+0*det([1,1,1;xy(nodes(j,:),:)'])*dt*f(sum(xy(nodes(j,:),:))/3,n*dt)/6;
     end
    %previous time step
    b=b+B*U(:,n-1);
    %Dirichlet Conditions
    u=sparse(size(xy,1),1);
    u(unique(Dirich)) = u_d(xy(unique(Dirich),:),n*dt);
    b=b-(dt*D*A+B-dt*C)*u;
    %computation of solution
    u(Freenodes)=(dt*D*A(Freenodes,Freenodes)+B(Freenodes,Freenodes)-dt*C(Freenodes,Freenodes))\b(Freenodes);
    U(:,n)=u;
  trisurf(nodes,xy(:,1),xy(:,2),U(:,n))
  xlabel('Horizontal Spatial Direction (meters)')
  ylabel('Vertical Spatial Direction (meters)')
  zlabel('Concentration of Contaminant (kg/m^2)')
  title(n)
  title(n)
  axis([0,10,0,10,0,1])
  
drawnow;
 end
 

