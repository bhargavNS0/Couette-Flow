function x = TDMAfunc(a,b,c,d,n)

%--calculation of P and Q---
P = zeros(1,n);
Q = zeros(1,n);
for i=1:n
    if i==1
        P(i) = -b(i)/a(i);
        Q(i) = d(i)/a(i);
    elseif i==n
        P(i) = 0;
        Q(i) = (d(i)-c(i)*Q(i-1))/(a(i)+c(i)*P(i-1));
    else
        P(i) = -b(i)/(a(i)+c(i)*P(i-1));
        Q(i) = (d(i)-c(i)*Q(i-1))/(a(i)+c(i)*P(i-1));
    end
end

%---Calculation of variables---
x = zeros(1,n);
for i=n:-1:1
    if i==n
        x(i) = Q(i);
    else
        x(i) = P(i)*x(i+1) + Q(i);
    end
end

end
