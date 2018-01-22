function [x, m , iter] = Chouliaras_assignment2_exercise2(f,df,x0,tol,nmax)
% G.C. Chouliaras
%This function determines the zero of a function using the Modified Newton's method.
%The inputs are: the function f, its derivative df, a starting point.
%Optional inputs are: the desired precision tol and the maximum number of
%iterations nmax. If they are not given as input, the default values are assigned.
%
%The outputs are: the approximated zero x, the estimate of the order of the
%zero m and the vector of all the iterations
%
%If the derivative of the function is zero, the function stops and a
%message informs the user.
%
%
%Example syntax:
%[x, m , iter] = Chouliaras_assignment2_exercise2(@(x) x^3 - 5, @(x) 3*x^2,5,0.00001,1500)

%check if tol exists, give default value
if ~exist('tol')
    tol=0.0001;
end

%check if nmax exists, give default value
if ~exist('nmax')
    nmax=1000;
end

%initial value x0
xnow = x0;

%initialize
niter = 0;
diff = tol+1; % to get into the loop

%set an initial value to M
M=1;

%evaluate the function and the derivative
fx = feval(f,xnow);
dfx = feval(df,xnow);


%create vector with approximated zeros
iter = zeros(1, nmax);

%set the first element equal to starting point x0
iter(1) = xnow;

% boolean variable check if we have more than 3 iterations
check = false;

%vector which contains the values of x
xvalues = zeros(1,3);

%the first value
xvalues(1) = xnow;

i = 2; %2 because we have the 1st (xnow)

%loop while maximum number of iterations is not reached, the difference
%between the last two values is larger than the desired precision, and the
%absolute value of the function evaluated at the last approximation is
%larger than the desired precision


while  niter <= nmax && diff >= tol && abs(fx) >= tol
    if check == false   %if we don't have 3 points
        niter = niter + 1;    % because niter =0 initially
        
        %check if derivative of x is 0
        if dfx == 0
            fprintf(['The derivative of f is 0']);
            break;
        end
        
        xnext = xnow - M * fx/dfx;  %use the modified newton's formula
        
        %difference between the last two values
        diff = abs(xnext - xnow);
        
        %update the points
        xnow = xnext;
        fx = feval(f,xnow);
        dfx = feval(df,xnow);
        
        xvalues(i) = xnow;
        i = i + 1;
        iter(niter + 1)= xnow;
        
        %check if we have 3 points
        
        if i > 3;
            check = true;
            i=2;
        end
        
        
        
    else
        
        %if we have 3 points, estimate m
        m = M * (xvalues(1) - xvalues(2))/ (xvalues(3) - 2* xvalues(2) + xvalues(1));
        
        %set M = m
        M = m;
        
        %update xvalues and set check to false, to use them again
        xvalues(1) = xnow;
        check = false;
        
        
    end
    
end

% if the maximum number of iterations is reached but the desired precision
% is not achieved, inform the user
if (niter==nmax && (abs(fx) > tol || diff > tol))
    fprintf(['The function stopped without converging to',...
        ' the desired tolerance because the maximum\n ',...
        'number of iterations was reached\n']);
end

%use the latest values to produce the output
x = xnow;
m = M;
iter = iter(1:niter + 1);



return