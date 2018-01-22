function sineapprox = Chouliaras_assignment1_exercise1(x)
% G.C. Chouliaras
%This function takes as an input a scalar or vector x and returns as output the sine of x
% as Taylor series expansion around 0 using a polynomial with 48 terms.
%
% Example:
%x = Chouliaras_assignment1_exercise1(0,4,10000)
%
sineapprox = zeros(size(x));
n=48; %determine the degree of the polynomial, number of iterations (48)

%%


for  z = 1:length(x) %loop through the length of vector x
    for i = 0:n %loop through the number of iterations n
      y(i + 1) = (-1)^i*x(z)^(2*i+1)/factorial(2*i+1); %calculate the terms
      % y
    end
    l = sum(y); % Sum all the increments
    sineapprox(z) = l; %Return the sum
end


    
%end 