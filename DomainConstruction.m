%% DomainConstruction.m Code
% This code reads in your data from the ChooseDomain.m file and builds your
% domain. You may run the CallAdvectionDiffusion.m file after running this
% code.
clear; clc;
% x=linspace(10,345,40);
% savex=length(x);
% y=linspace(10,495,40);
load mychosenpointsNeumann.mat
load mychosenpointsInterior.mat
load mychosenpointsDirichlet.mat
xy=[mychosenpointsInterior;mychosenpointsNeumann;mychosenpointsDirichlet];
x=xy(:,1);
y=xy(:,2);
savexsizemin = min(x);
saveysizemin = min(y);
savexsizemax = max(x);
saveysizemax = max(y);
save('savexsizemin.mat','savexsizemin')
save('saveysizemin.mat','saveysizemin')
save('savexsizemax.mat','savexsizemax')
save('saveysizemax.mat','saveysizemax')
%xy=[xy;mychosenpoints;mychosenpointsDirichlet];
% epsilon=10^-6;
% counteps=0;
%   for i=1:length(xy(:,1))
%     for j=1:length(xy(:,1))
%       
%         if abs(xy(i,1)-xy(j,1))<epsilon && abs(xy(i,2)-xy(j,2))<epsilon && i~=j
%             xy(j,:)=[0,0];
%             counteps=counteps+1;
%         end
%         
%     end
%   end
%  xysave=xy;
%  xy=xy(any(xy~=0,2),:); %Code from www.mathworks.com user Pekka Kumpulainen (March 12 2012).
 %xy=[0,0;xy];
 save('xy.mat','xy')
 nodes = delaunay(x,y); %elements3 is nodes
 save('nodes.mat','nodes')
 
 myvec=zeros(length(mychosenpointsNeumann(:,1)));
 for i=1:max(nodes(:))
     for j=1:length(mychosenpointsNeumann(:,1))
     
     if mychosenpointsNeumann(j,1) == xy(i,1) && mychosenpointsNeumann(j,2) == xy(i,2)
     myvec(i,1) = i;
     end
          
     end
 end
 myvec=myvec(myvec~=0);
 
 neumann=zeros(length(myvec)-1,2);
 for i=1:length(myvec)-1
   
   neumann(i,1)=myvec(i,1);   
   neumann(i,2)=myvec(i+1,1);
   
 end
 

 save('neumann.mat','neumann')
 
 
 myvec2=zeros(length(mychosenpointsDirichlet(:,1)));
 for i=1:max(nodes(:))
     for j=1:length(mychosenpointsDirichlet(:,1))
     
     if mychosenpointsDirichlet(j,1) == xy(i,1) && mychosenpointsDirichlet(j,2) == xy(i,2)
     myvec2(i,1) = i;
     end
          
     end
 end
 myvec2=myvec2(myvec2~=0);
 
 Dirichlet=zeros(length(myvec2)-1,2);
 for i=1:length(myvec2)-1
   
   Dirichlet(i,1)=myvec2(i,1);   
   Dirichlet(i,2)=myvec2(i+1,1);
   
 end
   Dirich=Dirichlet;
 Dirich=[neumann(end,end),Dirichlet(1,1);Dirich;Dirichlet(end,end),neumann(1,1)];
   
   
  save('Dirich.mat','Dirich')