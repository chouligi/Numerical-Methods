function determine_assignment3_exercise3

%this functions computes the convolution and the derivatives based on the
%FFT and also compares with the calculations in paper

f = @(x) sin(2*pi*x);
g = @(x) cos(2*pi*x);


N=2^4;
%I create linspace (0,1,N)
v = linspace(0,1,N+1);
%I dont get the last point because it is the same as the 1st one!
v(end)=[];

L = 1; %because I set 2*pi


f_fft = Chouliaras_assignment3_exercise3(f(v),'nor');
%calculate the coefficients

cf=f_fft./N;

g_fft = Chouliaras_assignment3_exercise3(g(v),'nor');

%coefficients for g
cg = g_fft./N;

%coefficients of convolution
conv_coef= L.*cf.*cg;


%derive the convolution
fft_convolution = Chouliaras_assignment3_exercise3(conv_coef*N,'inv');

%verify with calculated convolution in paper
real_conv = sin(pi*v).*cos(pi*v); 

N2 = floor(N/2);


%compute the coefficients
for index = 1:N2+1
    cf_derivative(index) = 2*pi/L *1i.* (index - 1) *cf (index);
end

for index = (N2+2):N
    cf_derivative(index) = 2*pi/L *1i.* (index - (N+1)) *cf (index);
end


%take the inverse to calculate the derivative
f_der = Chouliaras_assignment3_exercise3(cf_derivative,'inv')*N;

%compare with the derivative calculated in paper
real_deriv = 2*pi*cos(2 * pi*v);




%plot real vs approximated convolution


figure (1)

plot(real_conv,'r')
xlabel('x')
ylabel('convolution')
title('Convolution using FFT against the "exact" calculation')
hold on;
plot(real(fft_convolution),'--')
hold off;
legend(' "exact" ','using FFT')





%compute the coefficients for second derivative

for index = 1:N2+1
    cf_derivative2(index) = -(2*pi/L)^2* (index - 1)^2 *cf (index);
end

for index = (N2+2):N
    cf_derivative2(index) = -(2*pi/L)^2* (index - (N+1))^2 *cf (index);
end

%compute the second derivative by taking the inverse and multiplying by N
f_der2 = Chouliaras_assignment3_exercise3(cf_derivative2,'inv')*N;

%compare with the derivative computed in paper
real_deriv2 = -4*pi^2*sin(2*pi*v);



%plot real vs approximated derivative


figure (2)

plot(real_deriv,'r')
xlabel('x')
ylabel('derivative')
title('Derivative using FFT against the "exact" calculation')
hold on;
plot(real(f_der),'--')
hold off;
legend(' "exact" ','using FFT')


figure (3)

plot(real_deriv2,'r')
xlabel('x')
ylabel('2nd derivative')
title('2nd Derivative using FFT against the "exact" calculation')
hold on;
plot(real(f_der2),'--')
hold off;
legend(' "exact" ','using FFT')
