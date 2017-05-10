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

% %Second domain
load xysecondmesh.mat
load nodessecondmesh.mat
load neumannsecondmesh.mat
load Dirichsecondmesh.mat
load savexsizeminsecondmesh.mat
load saveysizeminsecondmesh.mat
load savexsizemaxsecondmesh.mat
load saveysizemaxsecondmesh.mat
load plotu3mpers.mat
load plotv3mpers.mat
load myxysecondmesh.mat

myvelocities=[reshape(plotu',length(plotu(:,1))^2,1),reshape(plotv',length(plotv(:,1))^2,1)];
nNodes=unique(nodes);
v=zeros(length(nNodes(:,1)),2);
epsilon=10;
for i=1:length(nNodes)
    for j=1:length(myxy(:,1))
        if abs(xy(i,1)-myxy(j,1))<epsilon && abs(xy(i,2)-myxy(j,2))<epsilon
            v(i,1)=myvelocities(j,1);
            v(i,2)=myvelocities(j,2);
        end
        
    end
end
v=v/10;%v/6000;
v=v;

%Initial Condition, velocity vector field, diffusivity term, maximum time
%to run, and time step
constant = (1/(sqrt(2*pi)^2*2));
f=@(x,y) constant*exp(-.0002*(x-200).^2-.0002*(y-400).^2);
%v=[-20;0];
%try high diffusivity with the neumann boundary condition pulling shit out
%on the left
%Diffusivity is when water is at 22 degrees celcius
D=1;%1.13e-3/24/60/60; %1;%http://www.isn.ucsd.edu/courses/beng221/problems/2012/BENG221_Project%20-%20Ao-Ieong%20Change%20Gu.pdf,D=10 and v=v*(1/2300) works with only point source Gaussian
T=800; dt=.1;

N=T/dt;

x=xy(:,1);
y=xy(:,2);
%elements3 data is nodes
Freenodes = unique(setdiff(nodes,Dirich)); %Freenodes
 A = sparse(size(xy,1),size(xy,1));
 B = sparse(size(xy,1),size(xy,1));
 C = sparse(size(xy,1),size(xy,1));
 %U = zeros(size(xy,1),N+1);
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
        mytempmat=xy(nodes(j,:),:);
        indices=[find(mytempmat(1,1)==xy(:,1) & mytempmat(1,2)==xy(:,2));find(mytempmat(2,1)==xy(:,1) & mytempmat(2,2)==xy(:,2));find(mytempmat(3,1)==xy(:,1) & mytempmat(3,2)==xy(:,2))];
        C(nodes(j,:),nodes(j,:))=C(nodes(j,:),nodes(j,:))+det([1,1,1;xy(nodes(j,:),:)'])*reshape([G*v(indices(1,1),:)';G*v(indices(2,1),:)';G*v(indices(3,1),:)'],3,3)*(1/6);       
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
          D*dt*g(sum(xy(neumann(j,:),:))/2,n*dt)/2;
       end
     
    %previous time step
    b=b+B*U(:,1);
    %Dirichlet Conditions
    u=sparse(size(xy,1),1);
    %u(unique(Dirich)) = u_d(xy(unique(Dirich),:),n*dt);
    u(unique(Dirich)) = u_dz(xy(unique(Dirich),:),n*dt);
    b=b-(dt*D*A+B-dt*C)*u;
    %computation of solution
    u(Freenodes)=(dt*D*A(Freenodes,Freenodes)+B(Freenodes,Freenodes)-dt*C(Freenodes,Freenodes))\b(Freenodes);
    U(:,2)=u;
    
    
  
          image=imread('domain.png');
        imshow(image)
         hold on
   trisurf(nodes,xy(:,1),xy(:,2),U(:,2))
   title(['t=',num2str(dt*n),'s'])
   xlabel('Horizontal Distance (meters)')
   ylabel('Vertical Distance (meters)')
   axis([savexsizemin,savexsizemax,saveysizemin,saveysizemax,0,1])
         view(2)
    colormap jet
  h = colorbar;
 ylabel(h, 'Concentration of Contaminant (kg per square meter)')
   drawnow;


   hold off
   U(:,1)=U(:,2);
   
 end

