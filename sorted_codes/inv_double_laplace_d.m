function  v_d = inv_double_laplace_d(N,T,n,K,r,x,G)
% 6 Compute f->v
% Discrete case
N = N;
T = T;
n = n;  %discretization of time
K = K;
r = r;
x = x;
delta = T/n;
P = expm(G*delta); %G->P
m1 = 10;
m2 = 15;
D = diag(x);

v_d = zeros(N,1);
A = 20;
t = K; %a form of strike price
rou = nthroot(exp(-A), 2*n);


% Inverse double laplace transform with Euler Approximation
%{
for j = 0:m1   %approximation
    sum3 = zeros(N,1);
    for jj = 0:m2+j    %summation in approximation
        theta = A/(2*t)-1i*pi/t-1i*(jj-1)*pi/t;
        sum2 = zeros(N,1);
        for k = -n:(n-1)    %no approximation
            z = rou*exp(1i*k*pi/n);
            f_d = (expm(theta*D) - z*P)\ones(N,1);
            L_d = 1/(theta^2)*f_d - 1/((theta^2)*(1-z)) + (x')/(theta*(1-z)*(1-z*exp(r*delta)));
            sum2 = sum2 + (-1)^k*L_d;
        end
        sum3 = sum3 + (-1)^jj*sum2;
    end
    v_d = v_d + factorial(m1)/(factorial(j)*factorial(m1-j))*2^(-m1)*sum3;
end
for j = 1:m1   %approximation
    sum3 = zeros(N,1);
    for jj = 0:m2+j    %summation in approximation
        theta = A/(2*t)-1i*pi/t-1i*(-jj+1)*pi/t;
        sum2 = zeros(N,1);
        for k = -n:(n-1)    %no approximation
            z = rou*exp(1i*k*pi/n);
            f_d = (exp(theta*D) - z*P)\ones(N,1);
            L_d = 1/(theta^2)*f_d - 1/((theta^2)*(1-z)) + (x')/(theta*(1-z)*(1-z*exp(r*delta)));
            sum2 = sum2 + (-1)^k*L_d;
        end
        sum3 = sum3 + (-1)^jj*sum2;
    end
    v_d = v_d + factorial(m1)/(factorial(j)*factorial(m1-j))*2^(-m1)*sum3;
end
v_d = exp(A/2)/(4*n*rou^n*t)*v_d;
%}
end
