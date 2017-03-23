%% AdvectionDiffusionFiniteElements.m Code
%This function solves the advection-diffusion equation on a specified two-dimensional domain utilizing the finite element method.
%The user can use the ChooseDomain.m and
%DomainConstruction.m files to build their domain.
%Diffusion code was drawn from
%https://www.math.hu-berlin.de/~cc/cc_homepage/download/1999-AJ_CC_FS-50_Lines_of_Matlab.pdf.
%Advection Code also inspired by the previously listed source but ultimately produced by
%Jordan Trinka of Carroll College in Helena, Montana. When running
%this code, expect numerical saw-tooth error for the first couple of trials
%due to chosen triangularization. The user will need to adjust their
%domain triangulation, initial conditions, and other parameters to get a
%a good numerical solution on a specified domain.
function [EndVec] = AdvectionDiffusionFiniteElementsConstant(xy,nodes,neumann,Dirich,savexsizemin,saveysizemin,savexsizemax,saveysizemax,f,v,D,T,dt)

N=T/dt;

x=xy(:,1);
y=xy(:,2);
%elements3 data is nodes
Freenodes = unique(setdiff(nodes,Dirich)); %Freenodes
 A = sparse(size(xy,1),size(xy,1));
 B = sparse(size(xy,1),size(xy,1));
 C = sparse(size(xy,1),size(xy,1));
% U = zeros(size(xy,1),N+1);
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
%Conservation of Contaminant
Conservation=zeros(N,1);
Conservation(1,1)=sum(U(:,1));
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
    b=b+B*U(:,n-1);
    %Dirichlet Conditions
    u=sparse(size(xy,1),1);
    u(unique(Dirich(7:end,:))) = u_d(xy(unique(Dirich(7:end,:)),:),n*dt);
    u(unique(Dirich(1:end-7,:))) = u_dz(xy(unique(Dirich(1:end-7,:)),:),n*dt);
    b=b-(dt*D*A+B-dt*C)*u;
    %computation of solution
    u(Freenodes)=(dt*D*A(Freenodes,Freenodes)+B(Freenodes,Freenodes)-dt*C(Freenodes,Freenodes))\b(Freenodes);
    U(:,n)=u;
    Conservation(n,1)=sum(U(:,n));
          image=imread('domain.png');
       imshow(image)
        hold on
  trisurf(nodes,xy(:,1),xy(:,2),U(:,n))
  title(dt*n)
  axis([savexsizemin,savexsizemax,saveysizemin,saveysizemax,0,1])
        view(2)
   colormap jet
  drawnow;
%   if n==430000 || n==200000 || n==100000 || n==50000
% frame = getframe(1);
%   im = frame2im(frame);
%   [Matrix,map] = rgb2ind(im,256);
%   Animationpic = sprintf('mathy_img_%d.png', n) ;
%  imwrite(Matrix,map, Animationpic, 'png');
%     
%   end

  hold off
 end
 EndVec=U(:,end); 
 
 figure
plot(Conservation)
end

