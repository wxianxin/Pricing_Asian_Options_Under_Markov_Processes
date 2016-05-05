function G = CTMC_MJD(N,r,d,x)

%Need constants sigma,r,d
sigma = .126349; 
alpha = -.390078; 
littlesigma = .338796;
lambda = .174814;
%MJD v(dy)
levy = lambda*(exp(2*alpha+littlesigma^2)*(exp(littlesigma^2)-1)+(exp(alpha+.5*littlesigma^2)-1)^2);

% 3 Compute Lambda Matrices
Lambda_J = zeros(N);
Lambda_D = zeros(N);

% 3.1 Compute Lambda_J matrix
for i=2:N-1
   for j=1:N
       if j~=i
          if j==1
          alphamax = (x(j+1)+x(j))/(2*x(i))-1;
          Lambda_J(i,j)=lambda*(normcdf(log(alphamax+1),alpha,littlesigma)-0);
          elseif j==N 
          alphamin = (x(j)+x(j-1))/(2*x(i))-1;  
          Lambda_J(i,j)=lambda*(1-normcdf(log(alphamin+1),alpha,littlesigma));
          else
          alphamin = (x(j)+x(j-1))/(2*x(i))-1;    
          alphamax = (x(j+1)+x(j))/(2*x(i))-1;
          Lambda_J(i,j)=lambda*(normcdf(log(alphamax+1),alpha,littlesigma)-normcdf(log(alphamin+1),alpha,littlesigma));
          end
       end
   end
end

for i=2:N-1
   Lambda_J(i,i) = -sum(Lambda_J(i,:));
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
      Jx=0;Jxx=0;
      for j=1:N
          Jx=Jx+Lambda_J(i,j)*(x(j)-x(i));
          Jxx=Jxx+Lambda_J(i,j)*(x(j)-x(i))^2;
      end
      B=zeros(3,1); 
      B(2)=(r-d)*x(i)-Jx;
      B(3)=x(i)^2*(sigma^2+levy)-Jxx;
      %solve the equation
      R =linsolve(A,B);
      Lambda_D(i,i-1)=R(1);
      Lambda_D(i,i)=R(2);
      Lambda_D(i,i+1)=R(3);
  end

% 4 Compute CTMC Markov process Generator G
G = Lambda_D + Lambda_J;
end