% Numerical analysis of Couette flow using Crank-Nicolson scheme in FDM

clear
clc

n = 20;                 % No. of elements
Re = 5000;              % reynolds number
plateVelocity = 1;      % in m/s
deltaY = 1/n;
y = 0:deltaY:1;
deltaT = 0.5;           % Choosen as per CFL condition for stability
u = zeros(1,n+1);       % inital condition
u(n+1) = plateVelocity; % boundary condition : plate velocity

m = input("Enter the nth timestep: ");
u0_5 = zeros(1,m);      % velocity at half height at each time step

%---------------------Numerical solution-------------------------------
for i=1:m 
    E = i*deltaT/(Re*deltaY^2);
    a = (1+E)*ones(1,n-1);
    b = -E/2*ones(1,n-1);
    c = -E/2*ones(1,n-1);
    d = zeros(1,n-1);
    for j=1:n-1
        d(j) = (1-E)*u(j+1) + E*(u(j+2)+u(j))/2;
    end
    d(n-1) = d(n-1) +E*u(n+1)/2;
    c(1) = 0;
    b(n-1) = 0;
    
    %----calculation of P and Q-------
    P = zeros(1,n-1);
    Q = zeros(1,n-1);
    
    for k=1:n-1
        if k==1
            P(k) = -b(k)/a(k);
            Q(k) = d(k)/a(k);
        elseif k==n-1
            P(k) = 0;
            Q(k) = (d(k)-c(k)*Q(k-1))/(a(k)+c(k)*P(k-1));
        else
            P(k) = -b(k)/(a(k)+c(k)*P(k-1));
            Q(k) = (d(k)-c(k)*Q(k-1))/(a(k)+c(k)*P(k-1));
        end
    end
    
    %---Calculation of variables---
    v = zeros(1,n-1);         % velocity from deltaY to n-1 deltaY
    for l=n-1:-1:1
        if l==n-1
            v(l) = Q(l);
        else
            v(l) = P(l)*v(l+1) + Q(l);
        end
    end

    u = zeros(1,n+1);         % Velocity at all nodes
    u(1) = 0;
    u(2:n) = v;
    u(n+1) = plateVelocity;
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

%------------------------Checking the steady state -----------------
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
% In this case, steady state achieved at 150th timestep



