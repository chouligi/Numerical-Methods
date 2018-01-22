function [x,v] = Chouliaras_assignment2_exercise1(f,p1,p2,N,tol)
% G.C. Chouliaras
%This function determines the zero of a function using the Secant method.
%The input is a function of a real variable f, two starting points p1 and
%p2,the maximum number of iteration steps N (optional), the desired
%precision tol (optional)
%The output is the calculated zero x and a vector of approximated values
%obtained in the iterations v.
%
%Example Syntax:
%[x,v] = Chouliaras_assignment2_exercise1(@(x) x^3 - 5,0,1,30,0.001)


%%
%p1 is p(n-2)
%p2 is p(n-1)

%make number of iteration optional
if ~exist('N')
    N=20;
end
if ~exist('tol')
    tol=10^(-4);
end

%vector to store the approximated values
v = [];

%the first approximation
p = p2 - (f(p2)*(p2 - p1))/(f(p2) - f(p1));

i = 2; %because we already have the firsts 2 values

%the first approximation is stored to vector v
v = [v p];

%loop while maximum number of iterations is reached or the desired
%precision is achieved
while (i <= N && abs(f(p)) > tol)
    
    %give the points their new valeus
    p1 = p2;
    p2 = p;
    
    %approximate again
    p = p2 - (f(p2)*(p2 - p1))/(f(p2) - f(p1));
    
    %store the new approximation to vector v
    v = [v p];
    
    %update i
    i = i + 1;
    
end
%if the while loop has ended without achieving the desired precision, a
%message informs the user
if (f(p) > tol)
    fprintf(['Secant stopped without converging to',...
        ' the desired tolerance because the maximum\n ',...
        'number of iterations was reached\n']);
end
%the approximated zero x is the last approximation p
x = p;



