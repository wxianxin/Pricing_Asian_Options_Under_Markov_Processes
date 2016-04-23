function  v_d = inv_laplace_d(N,T,n,K,r,x,G)
% 6 Compute f->v 
% Discrete case
N = N;
T = T;
n = n;
K = K;
r = r;
x = x;
delta = T/n;
P = expm(G*delta); %G->P
m1 = 10;
m2 = 15;
I = diag(ones(1,N));
D = diag(x);

v_d = zeros(N,1);
A = 20;
n = T/delta; %discretization of time
t = K*(n+1); %a form of strike price
rou = nthroot(exp(-A), 2*n);


theta = A/(2*t)-1i*pi/t-1i*0*pi/t;
v_d = 1/(theta^2)*((I/expm(D*theta))^(-1)*P)^n*(expm(theta*D)^(-1))*ones(N,1)-1/theta^2*ones(N,1)+x'/theta*(1-exp((n+1)*r*delta))/(1-exp(r*delta))