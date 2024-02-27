% Numerical analysis of Couette flow using Explicit scheme in FDM

clear
clc

n = 20;                 % No. of elements
Re = 5000;              % reynolds number
plateVelocity = 1;      % in m/s
deltaY = 1/n;
y = 0:deltaY:1;
deltaT = 6;           % Choosen as per CFL condition for stability
u = zeros(1,n+1);       % inital condition
u(n+1) = plateVelocity; % boundary condition : plate velocity

m = input("Enter the nth timestep: ");
u0_5 = zeros(1,m);      % velocity at half height at each time step
v = zeros(1,n+1);
v(n+1) = plateVelocity;

%---------------------Numerical solution-------------------------------
E = deltaT/(Re*deltaY^2);
for i=1:m 
    for j=2:n
        v(j) = E*u(j-1) + (1-2*E)*u(j) + E*u(j+1);
    end
    u = v;         % Velocity at all nodes
    u0_5(i) = u(n/2+1);       % if n is even, take (n+1)/2 if n is odd
    
end

%--------------------------Analytical solution-------------------
u1 = zeros(1,n+1);
for i=1:n+1
    u1(i) = y(i)*plateVelocity/y(n+1);
end

figure;
plot(u,y)
hold on
plot(u1,y,'r--')
xlabel("Velocity [m/s]");
ylabel("y [m]");
legend('Numerical','Analytical');
title("Velocity Distribution curve")
hold off

figure;
%t = deltaT:deltaT:m*deltaT;
plot(1:m,u0_5);
xlabel("\Delta t");
ylabel("Velocity [m/s]");

%------------------------Checking the steady -----------------
error = zeros(1,n+1);
sum = 0;
for i=2:n+1
    error(i) = u1(i)-u(i);
    if error(i) > 0.00001
        sum = sum + 1;
    end
end

if sum==0
    fprintf("Steady state achieved!"); 
else
    fprintf("Unsteady state prevails...")
end
% In this case, steady state achieved at 931th timestep