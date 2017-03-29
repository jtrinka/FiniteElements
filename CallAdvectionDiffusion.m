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
% load Dirichfirstmesh.mat
% load savexsizeminfirstmesh.mat
% load saveysizeminfirstmesh.mat
% load savexsizemaxfirstmesh.mat
% load saveysizemaxfirstmesh.mat
% load plotu3mpers.mat
% load plotv3mpers.mat
% load myxysecondmesh.mat

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
count=0;
for i=1:length(xy(:,1))
    
   if abs(xy(i,1)-150)<5
      count=count+1;
       wall(count,1)=i;
    
   end

end
% neumann1=neumann(1:122,:);
% neumann2=neumann(130:end,:);
% Dirich1=neumann(123:129,:);
% Dirich=[Dirich;Dirich1];
% neumann=[neumann1;neumann2];

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
WallPoint10Dec=AdvectionDiffusionFiniteElements(xy,nodes,neumann,Dirich,savexsizemin,saveysizemin,savexsizemax,saveysizemax,f,v,D-.1*D,T,dt,wall);


%% CallAdvectionDiffusion.m code
clearvars -except MyFirstEndVec; clc;
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
% load Dirichfirstmesh.mat
% load savexsizeminfirstmesh.mat
% load saveysizeminfirstmesh.mat
% load savexsizemaxfirstmesh.mat
% load saveysizemaxfirstmesh.mat
% load plotu3mpers.mat
% load plotv3mpers.mat
% load myxysecondmesh.mat

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

% neumann1=neumann(1:122,:);
% neumann2=neumann(130:end,:);
% Dirich1=neumann(123:129,:);
% Dirich=[Dirich;Dirich1];
% neumann=[neumann1;neumann2];

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

%Initial Condition, velocity vector field, diffusivity term, maximum time
%to run, and time step
constant = (1/(sqrt(2*pi)^2*2));
f=@(x,y) constant*exp(-.0002*(x-200).^2-.0002*(y-400).^2);
%v=[-20;0];
%try high diffusivity with the neumann boundary condition pulling shit out
%on the left
%Diffusivity is when water is at 22 degrees celcius
D=1;%1.13e-3/24/60/60; %1;%http://www.isn.ucsd.edu/courses/beng221/problems/2012/BENG221_Project%20-%20Ao-Ieong%20Change%20Gu.pdf,D=10 and v=v*(1/2300) works with only point source Gaussian
T=230; dt=.01;
MySecondEndVec=AdvectionDiffusionFiniteElements(xy,nodes,neumann,Dirich,savexsizemin,saveysizemin,savexsizemax,saveysizemax,f,v,D,T,dt);



%% CallAdvectionDiffusion.m code
clearvars -except MyFirstEndVec MySecondEndVec; clc;
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
% load Dirichfirstmesh.mat
% load savexsizeminfirstmesh.mat
% load saveysizeminfirstmesh.mat
% load savexsizemaxfirstmesh.mat
% load saveysizemaxfirstmesh.mat
% load plotu3mpers.mat
% load plotv3mpers.mat
% load myxysecondmesh.mat

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

% neumann1=neumann(1:122,:);
% neumann2=neumann(130:end,:);
% Dirich1=neumann(123:129,:);
% Dirich=[Dirich;Dirich1];
% neumann=[neumann1;neumann2];

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

%Initial Condition, velocity vector field, diffusivity term, maximum time
%to run, and time step
constant = (1/(sqrt(2*pi)^2*2));
f=@(x,y) constant*exp(-.0002*(x-200).^2-.0002*(y-400).^2);
%v=[-20;0];
%try high diffusivity with the neumann boundary condition pulling shit out
%on the left
%Diffusivity is when water is at 22 degrees celcius
D=1;%1.13e-3/24/60/60; %1;%http://www.isn.ucsd.edu/courses/beng221/problems/2012/BENG221_Project%20-%20Ao-Ieong%20Change%20Gu.pdf,D=10 and v=v*(1/2300) works with only point source Gaussian
T=230; dt=.001;
MyThirdEndVec=AdvectionDiffusionFiniteElements(xy,nodes,neumann,Dirich,savexsizemin,saveysizemin,savexsizemax,saveysizemax,f,v,D,T,dt);


%% CallAdvectionDiffusion.m code
clearvars -except MyFirstEndVec MySecondEndVec MyThirdEndVec; clc;
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
% load Dirichfirstmesh.mat
% load savexsizeminfirstmesh.mat
% load saveysizeminfirstmesh.mat
% load savexsizemaxfirstmesh.mat
% load saveysizemaxfirstmesh.mat
% load plotu3mpers.mat
% load plotv3mpers.mat
% load myxysecondmesh.mat

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

% neumann1=neumann(1:122,:);
% neumann2=neumann(130:end,:);
% Dirich1=neumann(123:129,:);
% Dirich=[Dirich;Dirich1];
% neumann=[neumann1;neumann2];

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

%Initial Condition, velocity vector field, diffusivity term, maximum time
%to run, and time step
constant = (1/(sqrt(2*pi)^2*2));
f=@(x,y) constant*exp(-.0002*(x-200).^2-.0002*(y-400).^2);
%v=[-20;0];
%try high diffusivity with the neumann boundary condition pulling shit out
%on the left
%Diffusivity is when water is at 22 degrees celcius
D=1;%1.13e-3/24/60/60; %1;%http://www.isn.ucsd.edu/courses/beng221/problems/2012/BENG221_Project%20-%20Ao-Ieong%20Change%20Gu.pdf,D=10 and v=v*(1/2300) works with only point source Gaussian
T=230; dt=.0001;
MyFourthEndVec=AdvectionDiffusionFiniteElements(xy,nodes,neumann,Dirich,savexsizemin,saveysizemin,savexsizemax,saveysizemax,f,v,D,T,dt);
save('MyFourthEndVec.mat','MyFourthEndVec')

%% Wall at 145 to 155
clear; clc;
load WallPoint10Inc
load WallPoint10Dec
load WallPointNoChange
plot(WallPointNoChange)
  hold on
  plot(WallPoint10Dec)
 hold on
plot(WallPoint10Inc)
title('Contaminant Flow through a Wall that Encloses 145 < x < 155')
xlabel('Time-Step (n)')
ylabel('Concentration of Contaminant (kg per square meter)')
legend D=1 D=0.9 D=1.1
 deltx=WallPoint10Inc(1900,1)-WallPointNoChange(1900,1);
 deltr=0.1;
 Sens=(deltx*1.1)/(WallPoint10Inc(1900,1)*deltr)
