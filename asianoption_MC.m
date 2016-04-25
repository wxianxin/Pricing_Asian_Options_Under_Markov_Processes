% Date: 20160424
% Arithmetic asian option

S0 =100;       % Price of underlying today
K = 100;       % Strike at expiry
sigma = 0.3;    % expected vol.
r = 0.0367;
T = 1;
m = 250; %Number of periods of calculation
dt = T/m;
N = 1000;
theta = 1;


S = zeros(1,m + 1);
S(1) = S0;
summation = 0;
for i = 1:N
    for ii = 2:(m+1)
        S(ii) = S(ii-1)*exp((r-0.5*sigma^2)*dt + sigma*sqrt(dt)*normrnd(0,1));
    end
    c = max(0,sum(S(2:m+1))/m - K);
    summation = summation + c;
end

C = summation/N

%{
    sigma_z = sigma*sqrt((2*m+1)/(6*(m+1)));
    rou = 0.5*((r-0.5*sigma^2)+sigma_z^2);
    
    d1 = (log(S0/K)+(rou+0.5*sigma_z^2))/(sqrt(T)*sigma_z);
    d2 = (log(S0/K)+(rou-0.5*sigma_z^2))/(sqrt(T)*sigma_z);
    c = exp(-r*T)*(S0*exp(rou*T)*normcdf(d1)-K*normcdf(d2));
%}
