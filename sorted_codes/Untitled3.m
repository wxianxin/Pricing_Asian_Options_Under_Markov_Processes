% Date: 20160424
% 1 Initialize
S_0 = 100;
N = 50; %Number of price states
T=1;    %time
n = 12; %Discretization of time
K = 100;    %strike price
r = 0.0367;
d = 0.00;
delta = T/n;

x = x_generator(S_0,N);
G = CTMC(N, r, d, x);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%555
A = 20;
A1 = 20;
A2 = 20;

t1 = K; %a form of strike price
t2 = T;%??
I = eye(N);
l1 = 3;
l2 = l1;
D = diag(x);
  %Truncation
v_dc = zeros(N,1);
for j = -nruns:nruns    
    theta = A1/(2*l1*t1)-1i*j*pi/l1/t1;
    sum2 = zeros(N,1);
    for k = -nruns:nruns
        mu = A2/(2*l2*t2)-1i*k*pi/l2/t2;
        f_dc = (theta*D + mu*I-G)\ones(N,1);
        L_c = 1/(theta^2)*f_dc - 1/((theta^2)*(mu)) + (x')/(theta*(mu)*(mu-r));
        sum2 = sum2 + exp(-1i*k*pi/l2)*L_c;
    end
    v_dc = v_dc +exp(-1i*j*pi/l1)*sum2;
end 
v_dc = exp(A/l1)/(4*l1*l2*t2*t1)*v_dc;
V_dc = exp(-r*T)/(T)*v_dc;