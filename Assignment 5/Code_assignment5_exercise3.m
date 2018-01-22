function [tvector,dvector,matrix] = Chouliaras_assignment5_exercise3(u0,L,T,S,N,s)
%G.C. Chouliaras
%
%This function solves the Korteweg-de Vries equation, which describes
%moderately high waves in a swallow channel.
%
%The input parameters are: a function handle u0 with the (almost) periodic initial
%conditions, the length of the interval L, the length of time T, the number
%of time steps S and optionally, the number of spatial discretization
%points N (default N = 100), and a parameter s which controls the speed of
%the animation (default s = 1/40).
%
%The output parameters are: a vector of times tvector, a vector of
%discretization points dvector and a matrix which contains the solution at
%those points in time and space. 
%Additionally, the function produces 2 pictures: a 3-D picture which shows
%the height of the wave u with respect to time and space and a 2-D
%animation which shows the solution in space through the passing of time.
%
%Example syntax:
%In order to solve the KDV equation with initial conditions u0(x) =
%f(x;lamda) = lamda/(2*(cosh(sqrt(lamda)*x/2)^2))for lamda =2 , length of interval 10,
%50 time steps, 100 spatial discretization points and length of time 20,
%one should type:
%f=@(x,lamda)  lamda./(2*(cosh(sqrt(lamda)*x/2).^2));
%
%[tvector,dvector,matrix] = Chouliaras_assignment5_exercise3(@(x) f(x,2),10,20,50,100)



%Check if input parameter for controlling animation speed exists, if not
%assign default value s = 1/40
if ~exist('s','var')
    s = 1/40;
end

%Check if there is optional input for the number of spatial discretization
%points N, if not assign default value N = 100
if ~exist('N','var')
    N = 100;
end

%take N discretization intervals, compute delta x
Dx = 2*L/N;
%discretize the spatial variable
dvector = -L:Dx:L;

%evaluate the initial condition
%we start from 2nd point because u0 = uN
init = feval(u0,dvector(2:end));

%Use sub-function korteweg to create a system of N ODE's to plug in the ode45 solver
[tvector,matrix] = ode45(@(t,init) korteweg(init,Dx,N),[0,T],init,odeset('AbsTol',1e-8,'RelTol',1e-8));


%Plotting%

%add one more column to the solution matrix in order to be able to plot
%with waterfall (vectors must have the same lengths, dvector has length 101)
matrix=[matrix(:,end),matrix];

%find the min and max elements of matrix
%nested min/max because min(matrix) returns the min element of each column,
%hence apply one more min to find the min of them
minel=min(min(matrix));
maxel=max(max(matrix));


%linspace generates S points between 1 and length(tvector)
%ceil rounds towards infinity


%we store the points we want to plot according to the time step S
points=ceil(linspace(1, length(tvector), S));
%the times corresponding to the points for plotting
timeplot = tvector(points);
%the solution at the requested points
matrixplot = matrix(points,:);


figure (2)

%3D Plot of u(x,t) in space and time
waterfall(dvector,timeplot,matrixplot)
xlabel('x')
ylabel('t')
zlabel('u')
title('3D plot of the evolution of the solution','fontsize',13);
set(gca,'fontsize',13);

figure(1)

%Animation showing the solution in space through the passing of time
for index=1:length(timeplot)
    plot(dvector,matrixplot(index,:))
    axis([-L,L,minel,maxel])
    xlabel('x')
    ylabel('u')
    title('2D plot of the solution','fontsize',13);
    set(gca,'fontsize',13);
    %draw repeatedly
    drawnow
    %control speed of animation with optional parameter s
    pause(s);
end






    function [f] = korteweg(u,dx,n)
        %This function uses central differences to compute the derivatives
        %in order to turn the PDE into a system of ODE's, and then it passes
        %this system into ode45 solver.
        
        %using centered differences compute the 1st derivatives of y
        %where y = (du/dx)^2 + 3u^2
        %initialize y
        y = zeros(n,1);
        
        %Case: i = 2,...,N-1
        index  = [2:1:n-1];
        y(index) = (u(index + 1) - 2*u(index) + u(index - 1))/(dx)^2 + 3*u(index).^2;
        
        %Case: i = 1
        y(1) = (u(2)- 2*u(1) +u(n))/(dx)^2 + 3*u(1)^2;
        
        %Case: i = N
        y(n) = (u(1) - 2*u(n)+u(n-1))/(dx)^2 + 3*u(n)^2;
        
        
        %Do the same for f = - dy/dx. This is the system we want to solve,
        %in order to solve the PDE.
        
        %initialize
        f = zeros(n,1);
        
        %Case: i = 2,...,N-1
        f(index) = - (y(index + 1) - y(index - 1))/(2*dx);
        
        %Case: i = 1
        f(1) = - (y(2) - y(n))/(2*dx);
        
        %Case: i = N
        f(n) = - (y(1) - y(n-1))/(2*dx);
    end
end
%test initial condition 1
%It shows traveling wave with constant velocity
%f=@(x,lamda)  lamda./(2*(cosh(sqrt(lamda)*x/2).^2));
%for lambda =2
%u0 = @(x) f(x,2)
%L = 10
%N = 100
%S = 150
%T = 20



%test initial condition 2, lamda1 = 14, lambda2 = 7
%Interaction between waves
%g = @(x) f(x - (5)/2, 14) + f(x + (5)/2,7);
%f=@(x,lamda)  lamda./(2*(cosh(sqrt(lamda)*x/2).^2));
%u0 = @(x) g(x)
%N = 100
%S = 500
%L = 5
%T = 4

%test initial condition 3, mu = 5, with T = 0.4 it shows breaking of waves,
% f = @(x) 5 * sin(pi*x/(100))
% T = 3
% S = 200
% N = 100
% L = 100





