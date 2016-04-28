% Date: 20160424
% 1 Initialize
S_0 = 1;
N = 50; %Number of price states
T = 1;	%time
n = 12; %Discretization of time
K = 1;    %strike price
r = 0.04;
d = 0.00;

x = x_generator(S_0,N);
G = CTMC(N, r, d, x);

v_c = qianyun_inv_double_laplace_c(N,T,n,K,r,x,G);
V_c = (exp(-r*T)/T) * v_c;

V_c(25)

%{
v_d = inv_double_laplace_d(N,T,n,K,r,x,G);
V_d = (exp(-r*T)/(n+1)) * v_d;
%}





