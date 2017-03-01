%% ChooseDomain.m Code
%This code allows a user to enter in any image of their choice which should
%be representitive of the domain they wish to solve the advection-diffusion
%equation on a two-dimensional domain by the finite element method. This is the first piece of code that should be ran when
%solving the advection-diffusion equation on a two-dimensional domain utilizing the finite element method. Then
%the user should run the DomainConstruction.m file and then finish by running the 
%CallAdvectionDiffusion.m file. When you run this code, your
%chosen image will appear and you first select which nodes will be Neumann
%boundary nodes, Dirichlet boundary nodes, and then free nodes. When
%outlining domain with boundary nodes, do not double over a node. Generate
%nodes in the sequential order you wish to have your boundaries formed
%between nodes. Also, do not select the same space twice with any node.
%This code will then save your node data in mychosenpointsNeumann.mat, mychosenpointsDirichlet.mat, and mychosenpointsInterior.mat.
%Run the DomainConstruction.m file after running this code completely.
%Code inspired by Dr. Eric Sullivan of Carroll College in Helena, MT.
%Neumann boundary
clear; clc;
A=imread('domain.png');
imshow(A)
clc;
fprintf('Click to choose Neumann nodes of domain. \n You will warned when the next chosen point will be a Dirichlet point. \n Then you will be warned again when the next chosen point will be Interior')
numberofpoints=200;
xNeumann = zeros(numberofpoints,1);
yNeumann = zeros(numberofpoints,1);
NeumannCounter=0;
for n=1:numberofpoints
    [x,y]= ginput(1);
   hold on
    plot(x,y,'r*')
    xNeumann(n,1)=x;
    yNeumann(n,1)=y;
    NeumannCounter=n
    
    if NeumannCounter+3==numberofpoints
       fprintf('\n Three Chosen Points until the next chosen point will be a Dirichlet Point')
       
    elseif NeumannCounter==numberofpoints
        fprintf('\nNext Chosen Point will be a Dirichlet Point')
    end
    
end
neumann=[xNeumann,yNeumann];
mychosenpointsNeumann=neumann;
%save as a .mat file
save('mychosenpointsNeumann.mat','mychosenpointsNeumann');

numberofpoints=10;
xDirichlet = zeros(numberofpoints,1);
yDirichlet = zeros(numberofpoints,1);
% DirichletCounter=0;
for n=1:numberofpoints
    [x,y]= ginput(1);
   hold on
    plot(x,y,'b*')
    xDirichlet(n,1)=x;
    yDirichlet(n,1)=y;
     DirichletCounter=n;
    
      if DirichletCounter+3==numberofpoints
       fprintf('\n Three Chosen Points until the next chosen point will be an Interior Point')
       
    elseif DirichletCounter==numberofpoints
        fprintf('\nNext Chosen Point will be an Interior Point')
      end

end
Dirichletpoints=[xDirichlet,yDirichlet];
mychosenpointsDirichlet=Dirichletpoints;
%save as a .mat file
save('mychosenpointsDirichlet.mat','mychosenpointsDirichlet');
%%
% A=imread('domain.png');
% imshow(A)
numberofpoints=900;
xInterior = zeros(numberofpoints,1);
yInterior = zeros(numberofpoints,1);
InteriorCounter=0;

for n=1:numberofpoints
    [x,y]= ginput(1);
   hold on
    plot(x,y,'k*')
    xInterior(n,1)=x;
    yInterior(n,1)=y;
    InteriorCounter=n
    
       
%     if InteriorCounter+3==numberofpoints
%        fprintf('\n Three More Chosen Points until Neumann')
%     elseif InteriorCounter==numberofpoints
%         fprintf('\n Next Chosen Point will be Neumann') 
%     end
end
Interior=[xInterior,yInterior];
mychosenpointsInterior=Interior;
%save as a .mat file
save('mychosenpointsInterior.mat','mychosenpointsInterior');

