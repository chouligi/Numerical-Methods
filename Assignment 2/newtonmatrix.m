function [zero,res,niter]=newtonmatrix(f,df,x0,tol,nmax)
%G.C.Chouliaras
%This function estimates the zero of a function using the Newton method.
%The inputs are the function f, its derivative f, the initial point x0, the
%desired precision tol and the maximum number of iterations nmax. The
%initial point x0 can be a scalar or a matrix.
%The outputs are the zero of the function, the residue and the number of
%iterations niter.
%
%Example Syntax:
%[zero,res,niter]=newtonmatrix(f,df,[2 3;3 4],0.0001,1000)


x = x0;
%evaluate the function and the derivative using arrayfun (x is matrix)
fx = arrayfun(f,x);
dfx = arrayfun(df,x);
niter = 0;
n = length(x0);
%make sure to get into while loop
diff = tol+ ones(n);
cond = 1;
%loop while x is not empty and number of iterations is smaller than the
%maximum
while ~isempty(x(cond)) && niter < nmax
    niter = niter + 1;
    diff = - fx./dfx;
    x = x  - (fx./dfx);
    diff = abs(diff);
    fx = arrayfun(f,x);
    dfx = arrayfun(df,x);
    cond = (diff >= tol);  % condition for checking if x is empty, x(-2 > 0) = []!
    %when it becomes FALSE then x(cond) = [] and algorithm stops
end
if niter > nmax
    fprintf(['Newton stopped without converging to',...
        ' the desired tolerance because the maximum\n ',...
        'number of iterations was reached\n']);
end
zero = x; res = fx;
return
