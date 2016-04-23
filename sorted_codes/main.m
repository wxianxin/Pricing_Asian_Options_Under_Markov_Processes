% Date: 20160422

% 1 Initialize
S_0 = 100;
N = 50; %Number of price states
T=1;    %time
n = 12; %Discretization of time
K = 100;    %strike price
delta = T/n;
r = 0.0367;
d = 0.00;

x = x_generator(S_0,N);
G = CTMC(N, r, d, x);

%v_d = inv_double_laplace_d(N,T,n,K,r,x,G);
v_d = ilaplace(ilaplace(L_d,j,t)k,n)
V_d = (exp(-r*T)/(n+1)) * v_d;

v_c = inv_double_laplace_c(N,T,n,K,r,x,G);
V_c = (exp(-r*T)/T) * v_c;



