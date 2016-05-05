% Date: 20160424
% 1 Initialize
S_0 = 1;
N = 1000; %Number of price states
T = 1;	%time
n = 250; %Discretization of time
K = 1;    %strike price
r = 0.04;
d = 0.00;
delta = T / n;

x = x_generator(S_0,N);
G = CTMC(N, r, d, x);

v_c = inv_double_laplace_c(N,T,n,K,r,x,G);
V_c = (exp(-r*T)/T) * v_c;

V_c(500)








