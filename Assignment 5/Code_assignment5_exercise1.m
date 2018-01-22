function [T,U] = Chouliaras_assignment5_exercise1(f, S, u0, theta, h)
%G.C. Chouliaras
%
%This function solves the initial value problem for a system of
%differential equations using the theta method.
%The theta method for the ODE u' = f(t,u) is given by:
% u(n+1) = u(n) + h[(1 - theta)*f(t(n),u(n)) + theta*f(t(n+1),u(n+1))]
%For theta > 0 the method is implicit, hence we rewrite the scheme as nonlinear problem and 
%we use multinewton to solve it.
%
%The inputs of the function are:
%a function handle f which denotes the rhs of the 
%differential equation, a vector S which includes the requested 
%points in time, a column vector u0 with the initial conditions,
%the parameter theta (between 0 and 1) and the step size h. If the vector S consists of more
%than two elements, the output is given for those points only. If S
%consists of 2 elements the function gives output for all calculated 
%points in time. 
%
%The output consists of column vector of times T and a matrix U in which
%every row has to the solution at the time corresponding to that row in T.
%
%Example syntax:
%
%In order to solve u' = -u, in [0 20], with initial condition u0=5, theta =
%0.5 (Crank Nicolson scheme) and step size h =0.01 one should type:
%
%[T,U] = Chouliaras_assignment5_exercise1(@(t,u) -u, [0 20], 5, 0.5, 0.01)
%f = @(t,u) -u;
%S = [0 20];
%u0 = 5;


%example
%f=@(t,x) cos(x)+t;
%S = [1,10];
%u0 = 1;
%[t,y]=ode45(f,S,u0);
%plot(t,y)


%create vector that includes the points in time that the
%solution must be calculated
tvector = S(1):h:S(end);

%create matrix s to include the solutions 
%s has a number of rows equal to the number of elements in u0 and 
%number of columns equal to the legth of tvector

%store the lengths of u0 and tvector
ulength = length(u0);
tlength = length(tvector);

%initialize the matrix s with zeros
s = zeros(ulength,tlength); 
%the first column corresponds to the given initial value
s(:,1)=u0;

%we have an implicit scheme so we have to rewrite the scheme and solve 
%G = 0

%define g
%u corresponds to the unknown u(n+1)
G = @(u) u - s(:,1) - h*((1 - theta)*f(tvector(1),s(:,1)) + theta*f(tvector(2),u));


%write G appropriately to input in numjac
Gnum= @(t,u) G(u);
DG = @(u) numjac(Gnum,0,u,Gnum(0,u),eps);

%loop from 1 to tlength - 1, at each step solve the nonlinear equation
%using Newton

for index = 1:tlength - 1
    G = @(u) u - s(:,index) - h*((1 - theta)*f(tvector(index),s(:,index)) + theta*f(tvector(index+1),u));
    
    %find the solution at the next point using newton, given as initial
    %value the previous solution
    
    %for Newton: tolerance = 0.0001 and number of iterations N =150   
    s(:,index+1) = multinewton(G,DG,s(:,index),0.0001,150);

    %if the length of S is larger than 2, we must provide only those points in time.
    %however these points are not precisely the instants where we computed
    %the solution, hence we interpolate using splines.
    
    if length(S) > 2
        %spline interpolation
        s_approx = spline(tvector,s,S);
        
        %give the values to the output 
        %make T a column vector
        T = transpose(S);
        U = s_approx';
        
    else
        %make T column vector
        T = transpose(tvector);
        U = s';
        
    end
    
end



