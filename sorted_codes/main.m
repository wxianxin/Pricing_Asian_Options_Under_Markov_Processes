% Date: 20160424
% 1 Initialize
S_0 = 100;
N = 50; %Number of price states
T = 1;	%time
n = 12; %Discretization of time
K = 100;    %strike price
r = 0.0367;
d = 0.00;
delta = T/n;
nruns = 200;

x = x_generator(S_0,N);
G = CTMC(N, r, d, x);

v_c = inv_double_laplace_c(N,T,n,K,r,x,G,nruns);
V_c = (exp(-r*T)/T) * v_c;


v_d = truncate_d(N,T,n,K,r,x,G);
V_dt100 = (exp(-r*T)/(n+1)) * v_d;

v_d = inv_double_laplace_d(N,T,n,K,r,x,G);
V_d = (exp(-r*T)/(n+1)) * v_d;






