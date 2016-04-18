% Date: 20160413

% 1 Initialize
S_0 = 100;
N = 1000;

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

%Need to define v(dy)
%Need code to compute integral y^2*v(dy) from -1 to infinity
%Need constants sigma,r,d

% 3 Compute Lambda Matrices
Lambda_J = zeros(N);
Lambda_D = zeros(N);

% 3.1 Compute Lambda_J matrix
for i = 2:(N-1)
    for j = 1:N
        if i ~= j
            if j == 1
                alphamin = 0;
                alphamax = rand*((x(j+1)/x(i))-(x(j)/x(i)))+(x(j)/x(i))-1;
            elseif j == N
                alphamin = rand*((x(j)/x(i))-(x(j-1)/x(i)))+(x(j-1)/x(i))-1;
                Lambda_J(i,j) = integral(fun,alphamin,Inf); %fun = v(dy)
            else
                alphamin = rand*((x(j)/x(i))-(x(j-1)/x(i)))+(x(j-1)/x(i))-1;
                alphamax = rand*((x(j+1)/x(i))-(x(j)/x(i)))+(x(j)/x(i))-1;
            Lambda_J(i,j) = integral(fun,alphamin,alphamax);
            end
        end
    end
    Lambda_J(i,i) = sum(Lambda_J(i,:));
end

% 3.2 Compute TRIDIAGONAL Lambda_D matrix
%Solve implicitly
Lambda_Jx = zeros(N);
Lambda_Jxx = zeros(N);
for i = 2:(N-1)
   A = zeros(3,N);
   for j = 1:N
       A(1,j) = 1;
       A(2,j) = x(j)-x(i);
       A(3,j) = (x(j)-x(i))^2;
   end
   for j = 1:N
   Lambda_Jx(i,j) = Lambda_J(i,j)*(x(j)-x(i));
   Lambda_Jxx(i,j) = Lambda_Jx(i,j)*(x(j)-x(i));
   end
    B = zeros(3,1);
    B(1,1)=0;
    %need levy constant (integral), sigma, r, d
    B(2,1)= (r-d)*x(i)-sum(Lambda_Jx(i,:));
    B(3,1)= x(i)*x(i)*(sigma(x(i))^2+levy)-sum(Lambda_Jxx(i,j));
    %solve
    R = linsolve(A,B);
    Lambda_D(i,:) = R;
end

% 4 Compute CTMC Markov process Generator G
G = Lambda_D + Lambda_J

% 5 Compute P
T=1;
n=N; %???????????
delta = T/N;
P = exp(delta*G)

% 6 Compute f
for c = 1:1000
    for cc = 1:1000 
        theta = 0.5 + 0.5i
        z = 0.5 + 0.5i
        mu = 0.5 + 0.5i
        D = diag(x);
        I = diag(ones(1,N));

        % 6.1 Discrete case
        f_d = inv(exp(1)^(theta*D) - z*P)

        % 6.2 Continuous case
        f_c = inv(theta*D + mu*I-G)

% 7.Compute L
r = 0.05
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
n = T/delta;
d = 0.00

% 9.1 Discrete case
V_d = (exp(-rT)/(n+1)) * v_d    %FUTURE CODE: Check the parameters

% 9.2 Continuous case
V_c = (exp(-rT)/T) * v_c    %FUTURE CODE: Check the parameters


