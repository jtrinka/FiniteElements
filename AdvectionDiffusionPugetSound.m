%% AdvectionDiffusionPugetSound.m Code
%This function runs the advection diffusion equation on a specified domain while utilizing the finite element method.
%The two-dimensional domain that we ran this code on originally was that of the portion of the Puget
%Sound that is close to Seattle. The user can use the ChooseDomain.m and
%DomainConstruction.m files to solve the advection-diffusion equation on a
%different, two-dimensional domain.
%Diffusion code was drawn from
%https://www.math.hu-berlin.de/~cc/cc_homepage/download/1999-AJ_CC_FS-50_Lines_of_Matlab.pdf.
%Advection Code also inspired by the previously listed source. When running
%this code, expect numerical saw-tooth error due to chosen triangularization. The user will need to adjust
%domain triangulation, initial conditions, and other parameters to get a
%a good numerical solution on a specified domain.
function [] = AdvectionDiffusionPugetSound(xy,nodes,neumann,Dirich,savexsizemin,saveysizemin,savexsizemax,saveysizemax,f,v,D,T,dt)

% load xy.mat
% load nodes.mat
% load neumann.mat
% load Dirich.mat
% load savexsizemin.mat
% load saveysizemin.mat
% load savexsizemax.mat
% load saveysizemax.mat
N=T/dt;

x=xy(:,1);
y=xy(:,2);
%elements3 data is nodes
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
%time steps
 for n=2:N+1
     b=sparse(size(xy,1),1);
%     %volume forces
     for j=1:size(nodes,1)
         b(nodes(j,:)) = b(nodes(j,:));%+0*det([1,1,1;xy(nodes(j,:),:)'])*dt*f(sum(xy(nodes(j,:),:))/3,n*dt)/6;
     end
     
     %Neumann conditions
       for j=1:size(neumann,1)
          b(neumann(j,:)) = b(neumann(j,:))+norm(xy(neumann(j,1),:)-...
          xy(neumann(j,2),:))*...
          dt*g(sum(xy(neumann(j,:),:))/2,n*dt)/2;
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
    image=imread('domain.png');
 imshow(image)
 hold on
  trisurf(nodes,xy(:,1),xy(:,2),U(:,n))
  title(n)
  axis([savexsizemin,savexsizemax,saveysizemin,saveysizemax,0,1])
  drawnow;
  hold off
 end
end
