% Date: 20160419

% 1 Initialize
S_0 = 100;
N = 10;
K = 101;    %strike price
T=1;    %tiem
r = 0.05;


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

n=N; %???????????
delta = T/N;
P = exp(delta*G); %????????????

% 6 Compute f->v
G = ones(N);    % FOR CODE BLOCK TEST, TO BE DELETED

D = diag(x);
V = zeros(N,1);
m1 = 10;
m2 = 15;

% 6.1 Discrete case
A = 20;
n = T/delta; %discretization of time
t = K/(n+1); %a form of strike price
rou = nthroot(exp(-A), 2*n);

for j = 0:m1   %approximation
    theta = A/(2*t)-1i*pi/t-1i*j*pi/t;
    sum3 = zeros(N,1);
    for jj = 0:m2+j    %summation in approximation
        sum2 = zeros(N,1);
        for k = -n:(n-1)    %no approximation
            z = rou*exp(1i*k*pi/n);
            f_d = inv(exp(theta*D) - z*P)*ones(N,1); % use a/b instead of a*inv(b)
            L_d = 1/(theta^2)*f_d - 1/((theta^2)*(1-z)) + (x')/(theta*(1-z)*(1-z*exp(r*delta_2)));
            sum2 = sum2 + (-1)^(k+1)*L_d;
        end
        sum3 = sum3 + (-1)^jj*sum2;
    end
    V = V + factorial(m1)/(factorial(j)*factorial(m1-j))*2^(-m1)*sum3;
end
for j = 0:m1   %approximation
    theta = A/(2*t)-1i*pi/t-1i*(-j)*pi/t;
    sum3 = zeros(N,1);
    for jj = 0:m2+j    %summation in approximation
        sum2 = zeros(N,1);
        for k = -n:(n-1)    %no approximation
            z = rou*exp(1i*k*pi/n);
            f_d = inv(exp(theta*D) - z*P)*ones(N,1); % use a/b instead of a*inv(b)
            L_d = 1/(theta^2)*f_d - 1/((theta^2)*(1-z)) + (x')/(theta*(1-z)*(1-z*exp(r*delta_2)));
            sum2 = sum2 + (-1)^(k+1)*L_d;
        end
        sum3 = sum3 + (-1)^jj*sum2;
    end
    V = V + factorial(m1)/(factorial(j)*factorial(m1-j))*2^(-m1)*sum3;
end
V = exp(A/2)/(4*n*rou^n*t)*V;

% 6.2 Continuous case
A1 = 20;
A2 = 20;
t1 = n;    %time
t2 = K/T;    %strike price
I = diag(ones(1,N));

for j = 0:m1    %approximation 2.1
    mu = A1/(2*t1)-1i*pi/t1-1i*j*pi/t1;
    sum3 = zeros(N,1);
    for jj = 0:m2+j    %summation in approximation 2.1
        sum2 = zeros(N,1);
        for jjj = 0:m1    %approximation 1.1
            sum1 = zeros(N,1);
            for jjjj = 0:m2+jjj    %summation in approximation 1.1
                theta = A2/(2*t2)-1i*pi/t2-1i*jjjj*pi/t2;
                f_c = inv(theta*D + mu*I-G)*ones(N,1); % use a/b instead of a*inv(b)
                L_c = 1/(theta^2)*f_c - 1/((theta^2)*mu) + (x')/(theta*mu*(mu-r));
                sum1 = sum1 + (-1)^jjjj*L_c;  
            end       
            sum2 = sum2 + factorial(m1)/(factorial(jjj)*factorial(m1-jjj))*2^(-m1)*sum1;
        end
        for jjj = 0:m1    %approximation 1.2
            sum1 = zeros(N,1);
            for jjjj = 0:m2+jjj    %summation in approximation 1.2
                theta = A2/(2*t2)-1i*pi/t2-1i*(-jjjj)*pi/t2;
                f_c = inv(theta*D + mu*I-G)*ones(N,1); % use a/b instead of a*inv(b)
                L_c = 1/(theta^2)*f_c - 1/((theta^2)*mu) + (x')/(theta*mu*(mu-r));
                sum1 = sum1 + (-1)^jjjj*L_c;  
            end       
            sum2 = sum2 + factorial(m1)/(factorial(jjj)*factorial(m1-jjj))*2^(-m1)*sum1;
        end
        sum3 = sum3 + (-1)^jj*sum2;
    end
    V = V + factorial(m1)/(factorial(j)*factorial(m1-j))*2^(-m1)*sum3;
end
for j = 0:m1    %approximation 2.2
    mu = A1/(2*t1)-1i*pi/t1-1i*(-j)*pi/t1;
    sum3 = zeros(N,1);
    for jj = 0:m2+j    %summation in approximation 2.2
        sum2 = zeros(N,1);
        for jjj = 0:m1    %approximation 1.1
            sum1 = zeros(N,1);
            for jjjj = 0:m2+jjj    %summation in approximation 1.1
                theta = A2/(2*t2)-1i*pi/t2-1i*jjjj*pi/t2;
                f_c = inv(theta*D + mu*I-G)*ones(N,1); % use a/b instead of a*inv(b)
                L_c = 1/(theta^2)*f_c - 1/((theta^2)*mu) + (x')/(theta*mu*(mu-r));
                sum1 = sum1 + (-1)^jjjj*L_c;  
            end       
            sum2 = sum2 + factorial(m1)/(factorial(jjj)*factorial(m1-jjj))*2^(-m1)*sum1;
        end
        for jjj = 0:m1    %approximation 1.2
            sum1 = zeros(N,1);
            for jjjj = 0:m2+jjj    %summation in approximation 1.2
                theta = A2/(2*t2)-1i*pi/t2-1i*(-jjjj)*pi/t2;
                f_c = inv(theta*D + mu*I-G)*ones(N,1); % use a/b instead of a*inv(b)
                L_c = 1/(theta^2)*f_c - 1/((theta^2)*mu) + (x')/(theta*mu*(mu-r));
                sum1 = sum1 + (-1)^jjjj*L_c;  
            end       
            sum2 = sum2 + factorial(m1)/(factorial(jjj)*factorial(m1-jjj))*2^(-m1)*sum1;
        end
        sum3 = sum3 + (-1)^jj*sum2;        
    end
    V = V + factorial(m1)/(factorial(j)*factorial(m1-j))*2^(-m1)*sum3;
end
V = exp((A1+A2)/2)/(4*t1*t2)*V;

% 7.Compute L
r = 0.05;
delta_2 = T/N;	%FUTURE CODE: What is the value of delta_2?

% 7.1 Discrete case
L_d = 1/(theta^2)*f_d - 1/((theta^2)*(1-z)) + (x')/(theta*(1-z)*(1-z*exp(r*delta_2)));

% 7.2 Continuous case
L_c = 1/(theta^2)*f_c - 1/((theta^2)*mu) + (x')/(theta*mu*(mu-r));

% 9 Compute option price
n = T/delta;
d = 0.00;

% 9.1 Discrete case
V_d = (exp(-rT)/(n+1)) * v_d    %FUTURE CODE: Check the parameters

% 9.2 Continuous case
V_c = (exp(-rT)/T) * v_c    %FUTURE CODE: Check the parameters


