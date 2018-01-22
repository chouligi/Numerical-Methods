function [T,posvel] = Chouliaras_assignment5_exercise2(m,x0,v0,Tin)
%G.C. Chouliaras
%
%This function shows the movement of planets of the solar system. Given a
%set of masses of the planets, their initial positions and velocities and
%the length of time over which the movement is computed, the function plots
%the orbits of the planets and also returns a vector of points in time and
%a matrix with the corresponding positions and velocities of the planets.
%
%The inputs of the function are:
%A vector m of length p with the masses of the planets, an nxp matrix of
%their initial positions x0, an nxp matrix of their initial velocities v0
%and the length of time Tin over which the movement of the planets is
%computed. For the inputs SI units are recommended in order to have
%appropriate results.
%
%The outputs of the function are:
%A vector of points in time T and a matrix posvel which includes the
%corresponding positions in the first half of the columns, and the corresponding
%velocities in the second half.
%
%Example syntax:
%
%In order to plot the Earth circling around the Sun for one year in n=2 dimensions one
%must do the following:
%
%The inputs:
%m =[1.9885e30;5.9724e24]; %kg
%x0 = [0, 1.496e11;0,0];%m
%v0 = [0,0;0,29800]; %m/s
%Tin = 365*24*60*60; %s
%
%Call the function:
%[T,posvel] = Chouliaras_assignment5_exercise2(m,x0,v0,Tin)


%store the dimensions (n) and the total number of planets (p)
%the dimensions can be n = 2 for a 2-D system, or n = 3 for a 3-D system.
[n,p] = size(x0);

%create a row vector of initial positions and velocities of the form (in case n = 2, p = 3):
%[x1 y1 x2 y2 x3 y3 ; vx1 vy1 vx2 vy2 vx3 vy3]'
%we make it in this form in order to approprietely insert it in ode45
initvector = [x0(:);v0(:)]';

%use sub-function f_handle as rhs for ODE and compute the positions and velocities of
%the planets
%we also add options with absolute and relative error tolerance in the ode45 for better performance

%[T, posvel] = ode45(@f_handle,[0 Tin],initvector,options,m,n,p);
[T, posvel] = ode45(@(t,initvector) f_handle(initvector,m,n,p),[0, Tin],initvector,odeset('AbsTol',1e-8,'RelTol',1e-8));

%the output of ode45 posvel includes the positions of the planets in the
%first n*p columns and the velocities in the rest of the columns.

%we store only the positions of the planets, in order to plot them
positions = posvel(:,1:(n*p));

%assign different color for each planet's orbit
colors = 'rbygkc';

%plotting
if n == 2
    
    for index = 1:p
        %if 2 dimensions
        %allocate space to store x,y positions in cells
        xpos = cell(p,1);
        ypos=  cell(p,1);
        
        %store the x and y positions
        xpos{index} = positions(:,2*index - 1);
        ypos{index} = positions(:,2*index);
        
        %plot the orbits
        figure(1)
        plot(xpos{index},ypos{index},'Color',colors(index));
        hold on
        plot(xpos{index}(1),ypos{index}(1),'yo','MarkerSize',20,'MarkerFaceColor',colors(index))
        set(gca,'fontsize',13);
        xlabel('x');
        ylabel('y');
    end
else
    for index = 1:p
        %if 3 dimensions, we include also a cell for the z-dimension
        xpos = cell(p,1);
        ypos = cell(p,1);
        zpos = cell(p,1);
        
        %store the x,y,z positions in the cells
        xpos{index} = positions(:,3*index - 2);
        ypos{index} = positions(:,3*index - 1);
        zpos{index} = positions(:,3*index);
  
        figure(1)
        
        plot3(positions(:,3*index - 2),positions(:,3*index - 1),positions(:,3*index),'Color',colors(index))
        axis equal
        hold on
      %plot initial conditions  
        plot3(xpos{index}(1),ypos{index}(1),zpos{index}(1),'yo','MarkerSize',20,'MarkerFaceColor',colors(index))
        set(gca,'fontsize',13);
        xlabel('x');
        ylabel('y');
        zlabel('z');
        
        
    end
   


end       

    function [output,A] = f_handle(thevector,m,n,p)
        %we use T as input, in order for ode45 to work
        
        %This sub-function takes as inputs the vector of initial positions
        %and velocities (thevector), the masses (m) , the number of
        %dimensions (n) and the number of planets (p). It returns the
        %vector output which consists the right-hand side for the ode45
        %solver.
        
        %the gravitational constant (in Nm^2 kg^(-2))
        G = 6.672*10^(-11);
        
        %store the positions from the vector of initial conditions in a row
        %vector pos
        pos = thevector(1:n*p);
        
        %reshape to a nxp matrix, whose elements are taken columnwise from row vector pos.
        %Each column corresponds to the position of a planet
        pos = reshape(pos,n,p);
        
        %create a nxp matrix with the total force on each planet in every column.
        
        %Use for loops to compute the total force on each planet and use
        %the total force to compute the accelerations
        
        %allocate space for the force and the acceleration
        F = zeros(n,p);
        Atemp = zeros(n,p);
        
        for index = 1:p
            for j  = 1:p
                %if the planets are different, compute the force between
                %planet index and planet j
                if index ~= j
                    %the force is computed as the size/magnitude of the force times
                    %the direction of the force.
                    Ftemp = (m(index) * m(j) * G *(pos(:,j) - pos(:,index))) / norm(pos(:,index) - pos(:,j),2)^3;
                else
                    %if index and j are the same, then the force is 0
                    Ftemp = 0;
                end
                %The total force applied on each planet is the sum of the forces
                %occur from the rest of the planets
                F(:,index) = F(:,index) + Ftemp;
            end
            %we compute accelarations using a = F/m and store them in a nxp matrix
            Atemp(:,index) = F(:,index)/m(index);
        end
        
        %turn acceleration matrix into a vector
        A = Atemp(:);
        
        %create system of differential equations.
        %we turn the system of nxp 2nd order ODE's to a 2*n*p 1st order system.
        
        %create a column vector output which indicates 2*n*p 1st order
        %ODE's
        %allocate space for the vector output
        output = zeros(2*n*p,1);
        
        for index = 1:n*p
            %the first n*p rows represent the velocity
            output(index) = thevector((n*p)+index);
            
            %the rest n*p rows represent the accelaration
            output((n*p) + index) = A(index);
            
        end
    end
end
%Test plots

%Sun - Earth 
%x0 = [0, 1.496e11;0,0]; %m
%v0 = [0,0;0,30000];%m/s
%m =[1.9891e30;5.97219e24]; %kg
%Tin = 365*24*60*60;

%Sun - Earth - Moon , n =2 dimensions 
%  m = [1.983537e30, 5.97219e24, 7.3477e22];
%x0 has form: [x1,x2,x3;y1,y2,y3]
%  x0=[0, 1.49e11, 1.49e11+4.05e8; 0,0,0];
%  v0 = [0,0,0; 0, 2.98e4, 2.98e4+964];
%  Tin = 365*24*60*60;


%sun-earth-saturn-mars
%m = [1.9891e30, 5.97219e24, 568e24,227.9e6]
%x0 = [0,149.6e9,1433.5e9,227.9e9;0,0,0,0;0,0,0,0]
%v0 = [0,0,0,0;0,29.8e3,9.7e3,24.1e3;0,0,0,0]
%Saturn needs 29.457 Earth years for a full circle around the sun
%Tin = 29.457*365*24*60*60;


%2 planets identical
%3D
%planets do not escape to infinity
%Tin = 365*24*60*60
%m = [1.9891e30,1.9891e30]
%x0=[1.496e+10,-1.496e+10;0,0;0,0];
%v0 = [0,0,;-29700,0;0,-29700];

%2D
%m = [1.9891e30,1.9891e30]
%x0=[1.496e+10,-1.496e+10;0,0]
%v0 = [0,0;29700,-29700]
%Tin = 365*24*60*60