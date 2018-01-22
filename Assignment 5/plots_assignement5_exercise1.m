function plots_assignement5_exercise1

%G.C. Chouliaras
%function for the plots of exercise 1. Each part must be executed
%separately.

%PLOT 1
%test function 1 for report
f = @(t,u) -u;
S = [0 20];
u0 = 5;
h=0.01;

%compare the methods (figure 1 in report)
figure(1)

[T,U] = Chouliaras_assignment5_exercise1(f, S, u0, 0.1, h);
plot(T,U,'b');
hold on
[T2,U2] = Chouliaras_assignment5_exercise1(f, S, u0, 1, h);
plot(T2,U2,'k');
hold on
[T3,U3] = Chouliaras_assignment5_exercise1(f, S, u0, 0.5, h);
plot(T3,U3,'g','LineWidth',7);
hold on
[t,y]=ode45(f,S,u0,odeset('AbsTol',1e-8,'RelTol',1e-8));
plot(t,y,'r');
legend('Forward Euler','Backward Euler','Crank-Nicolson','ode45');
        title('Solution of the ODE for 3 different thetas against ode45','fontsize',13)
        set(gca,'fontsize',13);
xlabel('t')
ylabel('u(t)')

hold off


%Plot 2
%test function 2 for report
%f inputs: number t, column vector u of length k, returns a column vector of length k, here k  = 2
f= @(t,u) [u(2);-u(1)];
S = [0 ,10*pi];
u0 = [1;0];
theta = 0.5;
h = 0.01;


[T,U] = Chouliaras_assignment5_exercise1(f, S, u0, theta, h);
[t,y]=ode45(f,S,u0,odeset('AbsTol',1e-8,'RelTol',1e-8));




%figure 2a
figure(2) 

plot(T,U,'*k');
hold on
plot(t,y,'g');
%legend('Crank-Nicolson','ode45');
legend('u1','u2','u1 ode45','u2 ode45')
title('Solution of system of ODEs with Crank-Nicolson and ode45','fontsize',13)
set(gca,'fontsize',13);
xlabel('t')
ylabel('u(t)')


figure(3)
[T,U] = Chouliaras_assignment5_exercise1(f, [30 330], u0, theta, 0.1);
[T2,U2] = Chouliaras_assignment5_exercise1(f, [30 330], u0, theta, 0.01);
[t,y]=ode45(f,[30 330],u0,odeset('AbsTol',1e-8,'RelTol',1e-8));

plot(T,U,'k');
hold on
plot(T2,U2,'r');
hold on
plot(t,y,'g');
legend('Crank-Nicolson','ode45');
legend('h=0.1','h=0.1','h=0.01','h=0.01','ode45','ode45')

%Plot 3
%Test function 3 for report
f = @(t,u) [2*u(2) + 8*u(1);-10*u(2) - 180*u(1)];
S = [0 5];
u0 = [40;-10];
theta = 0.5;
h = 0.01;

%the system of equations for 3a figure

[T,U] = Chouliaras_assignment5_exercise1(f, S, u0, theta, h);
[t,y]=ode45(f,S,u0,odeset('AbsTol',1e-8,'RelTol',1e-8));

figure(4)
plot(T,U,'*k');
hold on
plot(t,y,'g');
%legend('Crank-Nicolson','ode45');
legend('u1','u2','u1 ode45','u2 ode45')
title('Solution of system of ODEs with Crank-Nicolson and ode45','fontsize',13)
set(gca,'fontsize',13);
xlabel('t')
ylabel('u(t)')


