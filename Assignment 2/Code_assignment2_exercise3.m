function Chouliaras_assignment2_exercise3(V, discr, nsteps)
% G.C. Chouliaras
%The function applies Newton's method to the complex function f(z) = z^3 -
%1 and produces a picture of a rectangle in the complex plane, in which
%each point has a different color depending on which of the three zeros of
%the function the point is heading to.
%The input is: The vector of the two complex numbers V, a vector discr which
%indicates the number of discretization steps in both directions.
%Optional input is the number of iteration steps in Newton's method.
%There is no output.
%
%Example Syntax:
%Chouliaras_assignment2_exercise3([(-3 -3i), (3 + 3i)], [600, 600], 1000)




%assign default value to nsteps

if ~exist('nsteps')
    nsteps=1000;
end
%check the V to be appropriate vector! (to create appropriate rectangle)
if length(V) ~= 2 || length(discr) ~=2
    error('The length of both input vectors must be 2')
end

if V(1) == V(2)
    error('The input vectors must be different')
end

if imag(V(1)) ==0 && imag(V(2))==0
    error('The input vector must consist of two complex numbers')
end

%specify the function and its derivative
f = @(z) z^3 -1 ;
df = @(z) 3*z^2;

%assign the real and imaginary parts to variables
x1 = real(V(1));
y1 = imag(V(1));
x2 = real(V(2));
y2 = imag(V(2));



%calculate the appropriate step for x axis and y axis

xsteps = (abs(x2 - x1) / discr(1));
ysteps = (abs(y2 - y1) / discr(2));

%create vectors for starting values
%check values in order to create appropriate rectangle
if (x1 < x2)
    xvector = (x1:xsteps:x2);
else
    xvector = (x2:xsteps:x1);
end

if (y1 < y2)
    yvector = (y1:ysteps:y2);
else
    yvector = (y2:ysteps:y1);
end

%allocate space for the mesh/grid
z = zeros(length(xvector),length(yvector));

%create a vector with the roots
roots = [1 + 0i, -1/2 - sqrt(3)*1i/2, -1/2 + sqrt(3)*1i/2];




%replicate the xvector horizontally

xmatrix = repmat(xvector,length(xvector),1);

%replicate yvector vertically and multiply with i

ymatrix = repmat(yvector'.*1i,1,length(yvector));

%add the two matrices to create the grid
z = xmatrix + ymatrix;


%specify the precision for Newton's method
tol = 0.00001;

zero=newtonmatrix(f,df,z,tol,nsteps);



%the approximated points are not exactly equal to the roots, so we take the
%points that approximate the roots in the desired precision

%create 3 matrices with same size as matrix z, and include the roots
matrix1 = ones(length(xvector),length(yvector)).*roots(1);
matrix2 = ones(length(xvector),length(yvector)).*roots(2);
matrix3 = ones(length(xvector),length(yvector)).*roots(3);


%create 3 matrices with 1's and 0's, where 1 means that the approximated
%root is within the desired precision
red_matrix = abs(zero - matrix1) < tol ;
green_matrix = abs(zero - matrix2) < tol ;
blue_matrix = abs(zero - matrix3) < tol ;

%store the positions of the points that converge to each root
red_points = red_matrix == 1;
green_points = green_matrix == 1;
blue_points = blue_matrix == 1;

% plot the points
figure(1);


plot(z(red_points),'r.', 'MarkerSize',2);
hold on;
plot(z(green_points),'g.', 'MarkerSize',2);
hold on;
plot(z(blue_points),'b.', 'MarkerSize',2);
title('Convergence of Newton method to zeros of f(z) = z^3 - 1, depending on z0','fontsize',13);
xlabel('Re(z0)','fontsize',13);
ylabel('Im(z0)','fontsize',13);
set(gca,'fontsize',10);

hold off;

