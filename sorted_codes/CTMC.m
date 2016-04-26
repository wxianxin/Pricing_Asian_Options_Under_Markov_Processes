function G = CTMC(N,r,d,x)


%Need constants sigma,r,d
sigma = .126349; %???
alpha = -.390078; %mu is alpha in the paper
littlesigma = .338796;
lambda = .174814;
levy = lambda*((exp(2*alpha+(littlesigma^2))*(exp(littlesigma^2)-1))+(exp(alpha+.5*littlesigma^2)-1)^2);
%MJD v(dy)
fun = @(x,lambda,del,alpha) lambda./(x+1).*1./((2*3.1416)^.5)/del.*exp(-(log(x+1)-alpha).^2./2*(del.^2));
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
          Lambda_J(i,j)=lambda*(normcdf(alphamax,alpha,littlesigma)-normcdf(alphamin,alpha,littlesigma));
          elseif j==N 
          alphamin = rand*((x(j)/x(i))-(x(j-1)/x(i)))+(x(j-1)/x(i))-1;
          Lambda_J(i,j)=lambda*(1-normcdf(alphamin,alpha,littlesigma));
          else
          alphamin = rand*((x(j)/x(i))-(x(j-1)/x(i)))+(x(j-1)/x(i))-1;    
          alphamax = rand*((x(j+1)/x(i))-(x(j)/x(i)))+(x(j)/x(i))-1;
          Lambda_J(i,j)=lambda*(normcdf(alphamax,alpha,littlesigma)-normcdf(alphamin,alpha,littlesigma));
          end
       else
          Lambda_J(i,j) = -sum(Lambda_J(i,:));
       end
   end
end
%{
for i = 2:(N-1)
    for j = 1:N
        if i ~= j
            if j == 1
                alphamin = -1;
                alphamax = rand*((x(j+1)/x(i))-(x(j)/x(i)))+(x(j)/x(i))-1;
            elseif j == N
                alphamin = rand*((x(j)/x(i))-(x(j-1)/x(i)))+(x(j-1)/x(i))-1;
                Lambda_J(i,j) = integral(@(x)fun(x,.174814,.338796,-.390078),alphamin,Inf); %fun = v(dy)
            else
                alphamin = rand*((x(j)/x(i))-(x(j-1)/x(i)))+(x(j-1)/x(i))-1;
                alphamax = rand*((x(j+1)/x(i))-(x(j)/x(i)))+(x(j)/x(i))-1;
            Lambda_J(i,j) = integral(@(x)fun(x,.174814,.338796,-.390078),alphamin,alphamax);
            end
        end
    end
    Lambda_J(i,i) = -sum(Lambda_J(i,:));
end
%}

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
      levy = lambda*(exp(2*alpha+littlesigma^2)*(exp(littlesigma^2)-1)+(exp(alpha+.5*(littlesigma^2))-1)^2);
      B(3)=x(i)^2*(sigma^2+levy)-Jxx;
      %solve the equation
      R =linsolve(A,B);
      Lambda_D(i,i-1)=R(1);
      Lambda_D(i,i)=R(2);
      Lambda_D(i,i+1)=R(3);
  end
%{
%Solve implicitly
Lambda_Jx = zeros(N);
Lambda_Jxx = zeros(N);
for i = 2:(N-1)
   A = zeros(3,N);
   for j = 1:N
       A(1,j) = 1;
       A(2,j) = x(j)-x(i);
       A(3,j) = (x(j)-x(i))^2;
   end
   for j = 1:N
    Lambda_Jx(i,j) = Lambda_J(i,j)*(x(j)-x(i));
    Lambda_Jxx(i,j) = Lambda_Jx(i,j)*(x(j)-x(i));
   end
    B = zeros(3,1);
    B(1,1)=0;
    %need levy constant (integral), sigma, r, d
    B(2,1)= (r-d)*x(i)-sum(Lambda_Jx(i,:));
    %MJD sigma(xi)^2 = sigma^2
    B(3,1)= x(i)*x(i)*((sigma)^2+levy)-sum(Lambda_Jxx(i,:));
end
    %solve
R = linsolve(A,B);
for i = 2:(N-1)
    Lambda_D(i,:) = R;
end
%}

% 4 Compute CTMC Markov process Generator G
G = Lambda_D + Lambda_J;
end