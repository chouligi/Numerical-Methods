function [T,posvel,xpos,ypos] = plots_assignment5_exercise2(m,x0,v0,Tin)
%G.C. Chouliaras
%Function that creates the plot of the Earth around the sun.
%As inputs must be used the following:
%
%Sun - Earth
%x0 = [0, 1.496e11;0,0]; %m
%v0 = [0,0;0,29800];%m/s
%m =[1.9891e30;5.97219e24]; %kg
%Tin = 365*24*60*60; 

%store the dimensions (n) and the total number of planets (p)
[n,p] = size(x0);


%create vector of initial positions and velocities of form [x1 y1 x2 y2 x3
%y3 vx1 vy1 vx2 vy2 vx3 vy3]

initvector = [x0(:);v0(:)]';

%use f_handle as rhs for ODE and compute the positions and velocities of
%the planets
[T, posvel] = ode45(@(t,initvector) f_handle(initvector,m,n,p),[0, Tin],initvector,odeset('AbsTol',1e-8,'RelTol',1e-8));

%we store only the positions of the planets
positions = posvel(:,1:(n*p));


%assign different color for each planet's orbit
colors = 'rkgyc';

%plotting

for index = 1:p
    %if 2 dimensions
    if n == 2
        %allocate space to store x,y positions in cell
        xpos = cell(p,1);
        ypos=  cell(p,1);
        
        %store the x and y positions
        xpos{index} = positions(:,2*index - 1);
        ypos{index} = positions(:,2*index);
        

    else
        %if 3 dimensions
        xpos = cell(p,1);
        ypos = cell(p,1);
        zpos = cell(p,1);
        
        xpos{index} = positions(:,3*index - 2);
        ypos{index} = positions(:,3*index - 1);
        zpos{index} = positions(:,3*index);

        figure (1)
        plot3(xpos{index},ypos{index},zpos{index},'Color',colors(index));
        hold on
       plot3(xpos{index}(1),ypos{index}(1),zpos{index}(1),'yo','MarkerSize',20,'MarkerFaceColor',colors(index))

        axis equal
        %plot initial conditions  
        


    end

end 

%make plots for Earth - Sun System
%plot earth's initial position and final position

        figure (2)
        plot(xpos{2},ypos{2},'g')
        hold on
        plot(xpos{2}(1),0,'yo','MarkerSize',15,'MarkerFaceColor','b');
        hold on
        %plot the final position of Earth
        plot(xpos{2}(end),ypos{2}(end),'yo','MarkerSize',8,'MarkerFaceColor','r');
        hold on
        %plot the sun
        plot(0,0,'yo','MarkerSize',20,'MarkerFaceColor','y');
        legend('Orbit Earth','Initial position Earth','Final position Earth','Sun');
        title('Earth circling around the sun for 1 year','fontsize',13)
        set(gca,'fontsize',13);
        xlabel('x');
        ylabel('y');

 

        
      
