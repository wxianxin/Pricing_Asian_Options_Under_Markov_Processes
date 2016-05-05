% Date: 20160428
% Arithmetic asian option

S0 = 100;
K = 100;
sigma = 0.2;	
r = 0.04;
T = 1;
m = 1000;	%Number of periods of calculation
dt = T/m;
N = 1000;
theta = 1;


S = zeros(1,m + 1);
S(1) = S0;
summation = 0;
for i = 1:N
    for ii = 2:(m+1)
        S(ii) = exp(log(S(ii-1))+normrnd(0,sigma*sqrt(1/m)));
    end
    c = max(0,sum(S(2:m+1))/m - K);
    summation = summation + c;
end

C = exp(-r*T)*summation/N
