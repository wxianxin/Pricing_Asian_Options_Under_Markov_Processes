%Author: Steven Wang	Date: 20160413

% 1 Initialize
S_0 = 1;
N = 100;

% 2 Generate x
x(1) = 10^(-3) * S_0;
x(N) = 4 * S_0;
for k = 2:(N-1)
    if k <= N/2
        x(k) = S_0 + sinh((1-(k-1)/(N/2-1))*asinh(x(1)-S_0));
    else
        x(k) = S_0 + sinh(((k-N/2)/(N/2))*asinh(x(N)-S_0));
    end
end

% 3 Compute Lambda Matrices
Lambda_J = zeros(N);
Lambda_D = zeros(N);

% 3.1 Compute Lambda_J matrix
for i = 2:(N-1)
    for j = 1:N
        if i ~= j
            Lambda_J(i,j) = 0.1;    %FUTURE CODE: Compute the integral with the jump
        end
    end
    Lambda_J(i,i) = sum(Lambda_J(i,:));
end

% 3.2 Compute TRIDIAGONAL Lambda_D matrix

%FUTURE CODE

% 4 Compute CTMC Markov process Generator G
G = Lambda_D + Lambda_J

% 5 Compute P
T=1;
n=N; %???????????
delta = T/N;
P = exp(delta*G)

% 6 Compute f
theta = 0.5 + 0.5i    %FUTURE CODE: What is the value of theta?
z = 0.5 + 0.5i    %FUTURE CODE: What is the value of z?
mu = 0.5 + 0.5i    %FUTURE CODE: What is the value of mu?
D = diag(x);
I = diag(ones(1,N));

% 6.1 Discrete case
f_d = inv(exp(1)^(theta*D) - z*P)

% 6.2 Continuous case
f_c = inv(theta*D + mu*I-G)

% 7.Compute L
r = 0.02 %???????????
delta_2 = %FUTURE CODE: What is the value of delta_2?

% 7.1 Discrete case
L_d = 1/(theta^2)*f_d - 1/((theta^2)(1-z)) + x/(theta*(1-z)*(1-z*exp(r*delta_2)));

% 7.2 Continuous case
L_c = 1/(theta^2)*f_c - 1/((theta^2)*mu) + x/(theta*mu*(mu-r));

% 8 Compute inverse double transforms

% 8.1 Discrete case
v_d = %FUTURE CODE

% 8.2 Continuous case
v_c = %FUTURE CODE

% 9 Compute option price
n = T/delta; %???????
d = 0.02

% 9.1 Discrete case
V_d = (exp(-rT)/(n+1)) * v_d    %FUTURE CODE: Check the parameters

% 9.2 Continuous case
V_c = (exp(-rT)/T) * v_c    %FUTURE CODE: Check the parameters


