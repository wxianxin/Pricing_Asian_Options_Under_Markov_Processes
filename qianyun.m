clear all;
% 1 Initialize
S_0 = 100;
T=1;    %time
N = 50; %Number of price states
n = 12; %Discretization of time
K = 100;    %strike price
delta = T/n;
r = 0.0367;
d = 0.00;
mu = -.390078; %mu is alpha in the paper
littlesigma=.338796;
lambda=.174814;
sigma=.126349;

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

% 3 Compute Lambda Matrices
Lambda_J = zeros(N);
Lambda_D = zeros(N);
% 3.1 Compute Lambda_J matrix
for i=2:N-1
   for j=1:N
       if j~=i
          if j==1
          alphamin =0;
          alphamax = rand*((x(j+1)/x(i))-(x(j)/x(i)))+(x(j)/x(i))-1;
          Lambda_J(i,j)=lambda*(normcdf(alphamax,mu,littlesigma)-normcdf(alphamin,mu,littlesigma));
          elseif j==N 
          alphamin = rand*((x(j)/x(i))-(x(j-1)/x(i)))+(x(j-1)/x(i))-1;
          Lambda_J(i,j)=lambda*(1-normcdf(alphamin,mu,littlesigma));
          else
          alphamin = rand*((x(j)/x(i))-(x(j-1)/x(i)))+(x(j-1)/x(i))-1;    
          alphamax = rand*((x(j+1)/x(i))-(x(j)/x(i)))+(x(j)/x(i))-1;
          Lambda_J(i,j)=lambda*(normcdf(alphamax,mu,littlesigma)-normcdf(alphamin,mu,littlesigma));
          end
       else
          Lambda_J(i,j) = -sum(Lambda_J(i,:));
       end
   end
end

% 3.2 Compute TRIDIAGONAL Lambda_D matrix
  for i=2:N-1
      %MATRIX A
      A=ones(3);
      A(2,1)=x(i-1)-x(i);
      A(2,2)=x(i)-x(i);
      A(2,3)=x(i+1)-x(i);
      A(3,1)=(x(i-1)-x(i))^2;
      A(3,2)=(x(i)-x(i))^2;
      A(3,3)=(x(i+1)-x(i))^2;
      %MATRIX B
      B=zeros(3,1); 
      Jx=Lambda_J(i,i-1)*(x(i-1)-x(i))+Lambda_J(i,i)*(x(i)-x(i))+Lambda_J(i,i+1)*(x(i+1)-x(i));
      B(2)=(r-d)*x(i)-Jx;
      Jxx=Lambda_J(i,i-1)*(x(i-1)-x(i))^2+Lambda_J(i,i)*(x(i)-x(i))^2+Lambda_J(i,i+1)*(x(i+1)-x(i))^2;
      levy = lambda*(exp(2*mu+littlesigma^2)*(exp(littlesigma^2)-1)+(exp(mu+.5*(littlesigma^2))-1)^2);
      B(3)=x(i)^2*(sigma^2+levy)-Jxx;
      %solve the equation
      R =linsolve(A,B);
      Lambda_D(i,i-1)=R(1);
      Lambda_D(i,i)=R(2);
      Lambda_D(i,i+1)=R(3);
  end

% 4 Compute CTMC Markov process Generator G
G = Lambda_D + Lambda_J;

% 5 Compute P
P = expm(G*delta); %????????????