function thecurve = Chouliaras_assignment4_exercise1(f, df, x0, h0, N, minh, maxh)
%G.C. Chouliaras
%This function continues the solution of an equation f(x) = 0 where f: R^n
%-> R^(n-1). Given a function f, its derivative and a starting point, it
%calculates the solution curve of the function.
%The inputs are: the function f, the derivative df, a starting point x0 (must be a column vector), an
%initial step size h0, the maximum number of steps N and optionally the
%minimum step size minh and the maximum step size maxh. If the derivative of
%f is not given and the user inserts [] instead, then the function computes
%the derivative numerically using numjac.
%The output is a matrix of zeros of f (thecurve) which is a "vector" of
%points on the solution curve.
%
%
%Example syntax:
%In order to calculate the solution curve of f = @(x) x(1)^4 + x(2)^4 - 1,
%with df = @(x) [4*x(1)^3 4*x(2)^3] one should type:
%[thecurve] = Chouliaras_assignment4_exercise1(f, df, [1;0], 0.01, 5000, 0.0001, 2)


%check for minimum and maximum step size and provide values if they do not
%exist
if ~exist('minh','var')
    minh = 10^(-5);
end

if ~exist('maxh','var')
    maxh = 2;
end

if isempty(df) == 1
    %if df is [] we use numjac to compute the derivative numerically
    %in order to plug in the function to numjac we reform it
    fprime = @(t,x) f(x);
    df =@(x)  numjac(fprime, 0,x,f(x),eps);
end

%initialize
x_new = x0;

%find the tangent vector

v = null(df(x_new));

%take a small step in the direction of the tangent vector
xtilde = x_new + h0*v;
k = 0;

for i = 1:N
    
    %define function g which contains the two necessary conditions
    g = @(x)[f(x);v'*(x-xtilde)];
    %find the derivative of g
    dg = @(x) [df(x);v'];
    
    %give g as input to newton
    [x_newton,~,niter]=multinewton(g,dg,xtilde ,0.0001,1000);
    
    %Here the threshold must be very small in order to get good precision
    %If the number of iterations is smaller than 5 and f(x) very small,
    %we are at the correct direction so we accept the point and we increase
    %somewhat the step size
    if (niter < 5) & (f(x_newton) <= 10^(-10))
        
        h0 = min(1.1*h0, maxh);
        x_new = x_newton;
        k = k+1;
        thecurve(:,k) = x_new;
    else
        h0 = h0/2;
        %if the step size is smaller than its minimum value, the function
        %stops
        if h0 <minh
            error('Broke');
            break;
        end
    end
    
    
    vold = v;
    v = null(df(x_new));
    %check if the inner product is negative and if this is the case, change
    %sign of the tangent vector
    if dot(v,vold) < 0
        v = -v;
    end
    
    xtilde = (x_new + h0*v);
    
end


%Make the plot
% plot3(thecurve(1,:),thecurve(2,:),thecurve(3,:))
% title('Solution curve for f = [sin(x(1)^2) + log(x(2)); cos(x(3))]');
% set(gca,'fontsize',13);

end

