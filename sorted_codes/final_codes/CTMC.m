function G = CTMC(N,r,d,x)

N = N;
r = r;
d = d;
x = x;

%Need constants sigma,r,d
sigma = .7;
%alpha = -.390078; %mu is alpha in the paper
%littlesigma = .338796;

%lambda = .174814;

% 3 Compute Lambda Matrices

Lambda_D = zeros(N);

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
      B(2)=(r-d)*x(i);
      B(3)=x(i)^(1)*(sigma^2);
      %solve the equation
      R =linsolve(A,B);
      Lambda_D(i,i-1)=R(1);
      Lambda_D(i,i)=R(2);
      Lambda_D(i,i+1)=R(3);
  end


% 4 Compute CTMC Markov process Generator G
G = Lambda_D;
end