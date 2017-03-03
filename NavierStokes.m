%% Navier-Stokes from Cline
% at y = 0, yrange, all velocities = 0, due to boundaries
clear all
clc
clf
% load neumannfirstmesh.mat
% load xyfirstmesh.mat
load neumannsecondmesh.mat
load xysecondmesh.mat
neumann=unique(neumann);
mu = 10e-3;
xrange = 360; %meters, extent of domain in x direction
yrange = 520; %meters, extend of domain in y direction
xpoints = 100; %number of interior points in grid in x direction
ypoints = 100; %number of interior points in grid in y direction
tmax = 2; %seconds, duration of simulation
u0 = 3; %m/s, fan speed.
x = linspace(0,xrange,xpoints+2);
dx=x(2)-x(1);
y = linspace(0,yrange,ypoints+2);
dy=y(2)-y(1);
%The dt has to make sense
dt =3e-2;%sqrt(dx^2+dy^2)/340/2/30; % seconds; % seconds
tpoints = 120000;%ceil(tmax / dt) + 1;
[X,Y] = meshgrid(x,y);
checkx=X(1:end)';
checky=Y(1:end)';
myxy=[checkx,checky];
u = zeros(length(x),length(y),tpoints+1);
v = zeros(length(x),length(y),tpoints+1);
u_p = zeros(length(x),length(y));
v_p = zeros(length(x),length(y));
plusstuff=length(x)-82;

%initial conditions
u(85:95,70:80,:)=-u0;
u_p(85:95,70:80)=-u0;
v(85:95,70:80,:)=-u0;
v_p(85:95,70:80)=-u0;
u(1,:,:)=-u0;
u_p(1,:)=-u0;
v(1,:,:)=-u0;
v_p(1,:)=-u0;
 

rho0=1000;
epsx=2;
epsy=2;
for k = 1:tpoints
    for i = 2:length(x)-1
        for j = 2:length(y)-1
            
             
            %if ismember(i,2:64)==1 && j~=68 && ismember(i,71:82)==1
                
                dudt1(i,j) = -u(i,j,k)*(u(i+1,j,k) - u(i,j,k))/dx ...
                          -v(i,j,k)*(u(i,j+1,k) - u(i,j,k))/dy ...
                           + mu/rho0*((u(i+1,j,k)-2*u(i,j,k)+u(i-1,j,k))/dx^2 ...
                           +(u(i,j+1,k)-2*u(i,j,k)+u(i,j-1,k))/dy^2);
                       
            dvdt1(i,j) = -u(i,j,k)*(v(i+1,j,k) - v(i,j,k))/dx ...
                          -v(i,j,k)*(v(i,j+1,k) - v(i,j,k))/dy ...
                           + mu/rho0*((v(i+1,j,k)-2*v(i,j,k)+v(i-1,j,k))/dx^2 ...
                           +(v(i,j+1,k)-2*v(i,j,k)+v(i,j-1,k))/dy^2);
                       
                           
            u_p(i,j) = u(i,j,k) + dudt1(i,j) * dt;
            v_p(i,j) = v(i,j,k) + dvdt1(i,j) * dt;
            
                            
            
                     
           
            
            
            
        end
    end
%point diagonally up left
 %Working up
    u_p(85:95,70:80)=-u0;
    v_p(85:95,70:80)=-u0;
    
    u_p(85:97,70:75)=-u0;
    v_p(85:97,70:75)=-u0;
    
    u_p(85:99,60:69)=-u0;
    v_p(85:99,60:69)=-u0;
    %Working down
    u_p(85:93,80:87)=-u0;
    v_p(85:93,80:87)=-u0;
    
    u_p(70:84,88)=-u0;
    v_p(70:84,88)=-u0;
    
    u_p(70:82,89)=-u0;
    v_p(70:82,89)=-u0;
    
    u_p(70:80,90)=-u0;
    v_p(70:80,90)=-u0;
    
    u_p(70:78,91:92)=-u0;
    v_p(70:78,91:92)=-u0;
    
    u_p(70:77,93)=-u0;
    v_p(70:77,93)=-u0;
    
    u_p(61:73,94)=-u0;
    v_p(61:73,94)=-u0;
    
    u_p(55:61,95:99)=-u0;
    v_p(55:61,95:99)=-u0;
    
    u_p(55:61,100:101)=0;
    v_p(55:61,100:101)=-u0;
    
    u_p(48:55,95:100)=-u0;
    v_p(48:55,95:100)=-u0;
    
    u_p(43:47,95:98)=-u0;
    v_p(43:47,95:98)=-u0;
    
    u_p(40:42,95:98)=-u0;
    v_p(40:42,95:98)=-u0;
    
    u_p(38:40,95:99)=-u0;
    v_p(38:40,95:99)=-u0;
    
    u_p(37:39,100)=-u0;
    v_p(37:39,100)=-u0;
    
    u_p(33:38,101)=-u0;
    v_p(33:38,101)=-u0;
    
   
    
    
    
    for i = 2:length(x)-1
        for j = 2:length(y)-1
           
                
               
               dudt2(i,j) = -u_p(i,j)*(u_p(i,j) - u_p(i-1,j))/dx ...
                          -v_p(i,j)*(u_p(i,j) - u_p(i,j-1))/dy ...
                           + mu/rho0*((u_p(i+1,j)-2*u_p(i,j)+u_p(i-1,j))/dx^2 ...
                           +(u_p(i,j+1)-2*u_p(i,j)+u_p(i,j-1))/dy^2);
                      
          
            
            dvdt2(i,j) = -u_p(i,j)*(v_p(i,j) - v_p(i-1,j))/dx ...
                          -v_p(i,j)*(v_p(i,j) - v_p(i,j-1))/dy ...
                           + mu/rho0*((v_p(i+1,j)-2*v_p(i,j)+v_p(i-1,j))/dx^2 ...
                           +(v_p(i,j+1)-2*v_p(i,j)+v_p(i,j-1))/dy^2);
            
                           
            u(i,j,k+1) = u(i,j,k) + (dudt1(i,j)+dudt2(i,j)) * dt/2;
            v(i,j,k+1) = v(i,j,k) + (dvdt1(i,j)+dvdt2(i,j)) * dt/2;
                
          
            
            %Working up
u(85:95,70:80,k+1)=-u0;
v(85:95,70:80,k+1)=-u0;

u(85:97,70:75,k+1)=-u0;
v(85:97,70:75,k+1)=-u0;

u(85:99,60:69,k+1)=-u0;
v(85:99,60:69,k+1)=-u0;
%Working down
u(85:93,80:87,k+1)=-u0;
v(85:93,80:87,k+1)=-u0;

u(70:84,88,k+1)=-u0;
v(70:84,88,k+1)=-u0;

u(70:82,89,k+1)=-u0;
v(70:82,89,k+1)=-u0;

u(70:80,90,k+1)=-u0;
v(70:80,90,k+1)=-u0;

u(70:78,91:92,k+1)=-u0;
v(70:78,91:92,k+1)=-u0;

u(70:77,93,k+1)=-u0;
v(70:77,93,k+1)=-u0;

u(61:73,94,k+1)=-u0;
v(61:73,94,k+1)=-u0;

u(55:61,95:99,k+1)=-u0;
v(55:61,95:99,k+1)=-u0;

u(55:61,100:101,k+1)=0;
v(55:61,100:101,k+1)=-u0;

u(48:55,95:100,k+1)=-u0;
v(48:55,95:100,k+1)=-u0;

u(43:47,95:98,k+1)=-u0;
v(43:47,95:98,k+1)=-u0;

u(40:42,95:98,k+1)=-u0;
v(40:42,95:98,k+1)=-u0;

u(38:40,95:99,k+1)=-u0;
v(38:40,95:99,k+1)=-u0;
            
u(37:39,100,k+1)=-u0;
v(37:39,100,k+1)=-u0;

u(33:38,101,k+1)=-u0;
v(33:38,101,k+1)=-u0;




                  %for every i covers all possible j's, check it this way
                  
                  for counter = 1:length(neumann)
                      
                    
                      if abs(xy(neumann(counter,1),1)-X(j,i))<epsx && abs(xy(neumann(counter,1),2)-Y(j,i))<epsy
                      
                          u(i,j,k+1)=0;
                          v(i,j,k+1)=0;
                      
                      end
                     
                      
                  end

           
        end

    end


    plotu = u(:,:,k+1);
    plotv = v(:,:,k+1);
   image=imread('domain.png');
   imshow(image)
   hold on

    quiver(X,Y,plotu',plotv',1);
    title((k+1)*dt);
    drawnow;
    if k==10000 || k==20000 || k==30000 || k==40000 || k==50000 || k==60000 || k==70000 || k==80000 || k==100000
        
save('plotu.mat','plotu')
save('plotv.mat','plotv')
save('myxy.mat','myxy')
        
    end
       hold off
end
save('plotu.mat','plotu')
save('plotv.mat','plotv')
save('myxy.mat','myxy')