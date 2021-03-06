    %Uniform Domain Steady State Finite Element Method

    clear; clc;
    xmin = 0;
    xmax = 1;
    ymin = xmin;
    ymax = xmax;
    npt = 30;
    numintpt=npt;
    numintpts=npt;
    x = linspace(xmin,xmax,npt+2);
    dx = x(2)-x(1);
    y = linspace(ymin,ymax,npt+2);
    dy = y(2)-y(1);
    [x,y] = meshgrid(x,y);
    x=x';
    y=y';
    f = @(x,y) 20*exp(-((x-.5).^2+(y-.5).^2)/0.05);
    %Build A

    % Start by initializing the block structure of the stiffness matrix
    Z = zeros(numintpt,numintpt);
    for j=1:numintpt
        for k=1:numintpt
            C{j,k} = Z;
        end
    end
 
    % Now build the main diagonal block
    maindiag = 4*ones(numintpt,1);
    B1 = diag(maindiag);
    for j=2:numintpt
        B1(j-1,j)=-1;
        B1(j,j-1)=-1;
    end
 
    % Next build the upper and lower blocks
    B2 = diag(-ones(numintpt,1));
 
    % Fill the cell array with the appropriate blocks
    for j=1:numintpt
        C{j,j} = B1;
    end
    for j=2:numintpt
        C{j-1,j}=B2;
        C{j,j-1}=B2;
    end
 
    % build the matrix
    A = cell2mat(C);


    %create basis functions
    t = linspace(xmin,xmax,numintpt);
    w = linspace(ymin,ymax,numintpt);
    counter1=1; % I ADDED THIS IN
    counter2=1;
    for i=2:1:numintpt-1
        for j=2:1:numintpt-1
        
        %Region 1
        PQ1 = [x(i-1,j)-x(i,j);y(i-1,j)-y(i,j);-1];
        PR1 = [x(i,j-1)-x(i,j);y(i,j-1)-y(i,j);-1];
        Cross1 = cross(PQ1,PR1);
        R1 = @(t,w) ((((-Cross1(1,1)*(t-x(i,j))-Cross1(2,1)*(w-y(i,j)))/Cross1(3,1)))+1)*f(x(i,j),y(i,j));
        
        %Region 2
        PQ2 = [x(i-1,j)-x(i,j);y(i-1,j)-y(i,j);-1];
        PR2 = [x(i-1,j+1)-x(i,j);y(i-1,j+1)-y(i,j);-1];
        Cross2 = cross(PQ2,PR2);
        R2 = @(t,w) ((((-Cross2(1,1)*(t-x(i,j))-Cross2(2,1)*(w-y(i,j)))/Cross2(3,1)))+1)*f(x(i,j),y(i,j));
        
        %Region 3
        PQ3 = [x(i-1,j+1)-x(i,j);y(i-1,j+1)-y(i,j);-1];
        PR3 = [x(i,j+1)-x(i,j);y(i,j+1)-y(i,j);-1];
        Cross3 = cross(PQ3,PR3);
        R3 = @(t,w) ((((-Cross3(1,1)*(t-x(i,j))-Cross3(2,1)*(w-y(i,j)))/Cross3(3,1)))+1)*f(x(i,j),y(i,j));
        
        %Region 4
        PQ4 = [x(i,j+1)-x(i,j);y(i,j+1)-y(i,j);-1];
        PR4 = [x(i+1,j)-x(i,j);y(i+1,j)-y(i,j);-1];
        Cross4 = cross(PQ4,PR4);
        R4 = @(t,w) ((((-Cross4(1,1)*(t-x(i,j))-Cross4(2,1)*(w-y(i,j)))/Cross4(3,1)))+1)*f(x(i,j),y(i,j));
        
        %Region 5
        PQ5 = [x(i+1,j)-x(i,j);y(i+1,j)-y(i,j);-1];
        PR5 = [x(i+1,j-1)-x(i,j);y(i+1,j-1)-y(i,j);-1];
        Cross5 = cross(PQ5,PR5);
        R5 = @(t,w) ((((-Cross5(1,1)*(t-x(i,j))-Cross5(2,1)*(w-y(i,j)))/Cross5(3,1)))+1)*f(x(i,j),y(i,j));
        
        %Region 6
        PQ6 = [x(i+1,j-1)-x(i,j);y(i+1,j-1)-y(i,j);-1];
        PR6 = [x(i,j-1)-x(i,j);y(i,j-1)-y(i,j);-1];
        Cross6 = cross(PQ6,PR6);
        R6 = @(t,w) ((((-Cross6(1,1)*(t-x(i,j))-Cross6(2,1)*(w-y(i,j)))/Cross6(3,1)))+1)*f(x(i,j),y(i,j));
        
        %Region 1 line
        m1 = (y(i,j-1)-y(i-1,j))/(x(i,j-1)-x(i-1,j));
        b1 = (y(i-1,j)-m1*(x(i-1,j)));
        R1line = @(t) m1*t+b1;
        
        %Region 2/3 line
        m23 = (y(i,j)-y(i-1,j+1))/(x(i,j)-x(i-1,j+1));
        b23 = (y(i-1,j+1)-m23*(x(i-1,j+1)));
        R23line = @(t) m23*t+b23;
        
        
        %Region 4 line
        m4 = (y(i+1,j)-y(i,j+1))/(x(i+1,j)-x(i,j+1));
        b4 = (y(i,j+1)-m4*(x(i,j+1)));
        R4line = @(t) m4*t+b4;
        
        %Region 5/6 line
        m56 = (y(i+1,j-1)-y(i,j))/(x(i+1,j-1)-x(i,j));
        b56 = (y(i,j)-m56*(x(i,j)));
        R56line = @(t) m56*t+b56;
        
        %Region 1 integral
        R1int = integral2(R1,x(i-1,j),x(i,j),R1line,y(i,j));
        
        %Region 2 integral
        R2int = integral2(R2,x(i-1,j),x(i,j),y(i,j),R23line);
        
        %Region 3 integral
        R3int = integral2(R3,x(i-1,j),x(i,j),R23line,y(i,j+1));
        
        %Region 4 integral
        R4int = integral2(R4,x(i,j),x(i+1,j),y(i,j),R4line);
        
        %Region 5 integral
        R5int = integral2(R5,x(i,j),x(i+1,j),R56line,y(i,j));
        
        %Region 6 integral
        R6int = integral2(R6,x(i,j),x(i+1,j),y(i,j-1),R56line);
      
        
            H=R1int+R2int+R3int+R4int+R5int+R6int;
        
    RHS(i+1,j+1)=H;
            counter1 = counter1+1;
             counter2 = counter2+1;
        end
    end

    RHS = reshape(RHS',npt^2,1);

    uint = A\RHS;

    uint = reshape(uint,npt,npt);


    bound1 = zeros(numintpt,1);
    utot1 = zeros(numintpt,numintpt+2);
    utot1(:,(2:end-1)) = uint;
    bound2 = zeros(1,numintpt+2);
    utot2 = [bound2;utot1;bound2];

    surf(x,y,utot2)

