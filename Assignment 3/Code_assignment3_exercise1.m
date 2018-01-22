function [l_best,c_best,residue] = Chouliaras_assignment3_exercise1(data, g)
% G.C. Chouliaras
%This function fits a linear combination of exponential function to data.
%The inputs are: data is a m x 2 matrix and g is the vector of initial
%guesses.
%
%The outputs are: l_best is the vector of optimal expoments lamda_i, c_best
%is the vector of optimal constants C_i,residue is the residual

%
%
%Example syntax: [l_best,c_best,residue] = Chouliaras_assignment3_exercise1(data1, [0.1 0.2 0.3])
%where data1 is the first data set of expo-examples.mat

%number of exponents follows from vector of guesses
n = length(g);

%allocate space for vectors with points x's and y's
xdata = zeros(n,1);
ydata = zeros(n,1);

xdata = data(:,1);
ydata = data(:,2);


%create anonymous function A in which the matrix A is constructed
A  = @(l) (exp(1).^(xdata*l));

%calculate ci's, length of C is the same as length of the guesses
C = zeros(n,1);

%make coefficients C function of lamdas! l corresponds to g
C =  @(l) pinv(A(l))*ydata;

%function f, depends on lamda, f1,f2,...fn
f = @(l) A(l) * C(l);

%compute the residue, but take the squared norm
res = @(l) norm(f(l) - ydata)^2;

%find optimal lambdas by minimizing the residue
l_best = fminsearch(res,g);

%find optimal ci's using the optimal lambdas!
c_best = pinv(A(l_best))*ydata;

%the residue
residue = norm(A(l_best)*c_best - ydata)^2;


fit=zeros(length(xdata),1);

%create the fitted function

fit = A(l_best)*c_best;

%plot the fitted function against the data
figure (1);
plot(xdata,ydata,'r');
hold on;
plot(xdata,fit,'k.')
hold off;


end
