function v_d = truncate_d(N,T,n,K,r,x,G)
NN = 100;
N = N;
T = T;
n = n;
K = K;
r = r;
x = x;
delta = T/n;
P = expm(G*delta); %G->P
D = diag(x);

v_d = zeros(N,1);
A = 20;
n = T/delta; %discretization of time
t = K; %a form of strike price
rou = nthroot(exp(-A), 2*n);


v_d = 0;
for j = 1:NN
    sum1 = 0;
    theta = A/(2*t)-1i*pi/t-1i*(-j+1)*pi/t;
    for k = -n:(n-1)    %no approximation
        z = rou*exp(1i*k*pi/n);
        f_d = (exp(theta*D) - z*P)\ones(N,1);
        L_d = 1/(theta^2)*f_d - 1/((theta^2)*(1-z)) + (x')/(theta*(1-z)*(1-z*exp(r*delta)));
        sum1 = sum1 + (-1)^k*L_d;
    end
    v_d = v_d + sum1;

end