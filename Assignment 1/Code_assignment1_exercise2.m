function [t,x] = Chouliaras_assignment1_exercise2(m,b,a,x0,v0)

% G.C. Chouliaras
%This function solves m*x'' + b*x' +a*x = 0 with initial conditions x(0) =
%x0 and x'(0)=v0. The inputs are the parameters: m > 0, b >= 0 ,a>= 0 and the initial
%conditions x0, v0 (real numbers, can be scalars or vectors of the same size). 
%The output is a vector t (time) and a matrix x (solution).
%The size of x depends on the size of the vectors x0, v0
%
%If m,b,a are not appropriate an error is returned.
% Example:
%[t,x] = Chouliaras_assignment1_exercise2(2,1,3,[4 3],[1 3])

%%
%determine the interval [0,T] for which the ode is solved
tstart = 0;
tend = 50;
step = 0.0001;

t = (tstart:step:tend)' ;
n=length(t);

d1 = length(x0);
x = zeros(d1,n); %create matrix in which the solutions will be stored


if (m > 0 && b >= 0 && a >= 0) %if the input is appropriate
    
    D = b^2 - 4*m*a; %calculate determinant
    
    if D > 0   %case D > 0
        r1 = (-b - sqrt(D))/(2*m); %calculate roots
        r2 = (-b + sqrt(D))/(2*m);
        
        A = x0 - (v0 - x0.*r1)/(r2 - r1); %calculate constants A and B (if x0, v0 are vectors, these are vectors too)
        B = (v0 - x0.*r1)/(r2 - r1);
        
        x = A'*exp(t'.*r1) + B'*exp(t'.*r2);%calculate the solution

    elseif D == 0 % case D=0
        r = -b/(2*m);
        
        A = x0;
        B = v0 - x0.*r;
       
       
            x = A'*exp(t'.*r) + (t.*exp(t.*r)*B)'; 
         
    else
        
        l = -b/(2*m);
        omega = sqrt(4*m*a - b^2)/(2*m);
        
        A = x0;
        B = (v0-x0.*l)/omega;
        
      
            x = A'*(exp(t.*l).*cos(t.*omega))' + B'*(exp(t.*l).*sin(t.*omega))';
       
        
    end
    
else
    % eror message in case the solution is not appropriate
    error('One of these mistakes in input: m <= 0 or b < 0 or a < 0') 
    
end