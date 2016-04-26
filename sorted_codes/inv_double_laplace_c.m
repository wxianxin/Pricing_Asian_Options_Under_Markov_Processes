function  v_c = inv_double_laplace_c(N,T,n,K,r,x,G,nruns)
% 6 Compute f->v
% Continuous case
m1 = 10;
m2 = 15;
I = eye(N);
D = diag(x);
v_c = zeros(N,1);
A1 = 20;
A2 = 20;
A = 20;
t1 = K;    %strike price
t2 = T;    %time
l1 = 3;
l2 = l1;

%Truncation
for j = -nruns:nruns    
    theta = A/(2*l1*t1)-1i*j*pi/l1/t1;
    sum2 = zeros(N,1);
    for k = -nruns:nruns
        mu = A2/(2*l2*t2)-1i*k*pi/l2/t2;
        f_dc = (theta*D + mu*I-G)\ones(N,1);
        L_c = 1/(theta^2)*f_dc - 1/((theta^2)*(mu)) + (x')/(theta*(mu)*(mu-r));
        sum2 = sum2 + exp(-1i*k*pi/l2)*L_c;
    end
    v_c = v_c + exp(-1i*j*pi/l1)*sum2;
end
v_c = exp(A/l1)/(4*l1*l2*t2*t1)*v_c;

%v_c = exp(0.5*A1/l1)/(2*l1*t1) * v_c;
%Double inverse Laplace transform with Euler Approximation
%{
for j = 0:m1    %approximation 2.1
    sum3 = zeros(N,1);
    for jj = 0:m2+j    %summation in approximation 2.1
        mu = A1/(2*t1)-1i*pi/t1-1i*(jj-1)*pi/t1;
        sum2 = zeros(N,1);
        for jjj = 0:m1    %approximation 1.1
            sum1 = zeros(N,1);
            for jjjj = 0:m2+jjj    %summation in approximation 1.1
                theta = A2/(2*t2)-1i*pi/t2-1i*(jjjj-1)*pi/t2;
                f_c = (theta*D + mu*I-G)\ones(N,1);
                L_c = 1/(theta^2)*f_c - 1/((theta^2)*mu) + (x')/(theta*mu*(mu-r));
                sum1 = sum1 + (-1)^jjjj*L_c;  
            end       
            sum2 = sum2 + factorial(m1)/(factorial(jjj)*factorial(m1-jjj))*2^(-m1)*sum1;
        end
        for jjj = 1:m1    %approximation 1.2
            sum1 = zeros(N,1);
            for jjjj = 0:m2+jjj    %summation in approximation 1.2
                theta = A2/(2*t2)-1i*pi/t2-1i*(-jjjj+1)*pi/t2;
                f_c = (theta*D + mu*I-G)\ones(N,1);
                L_c = 1/(theta^2)*f_c - 1/((theta^2)*mu) + (x')/(theta*mu*(mu-r));
                sum1 = sum1 + (-1)^jjjj*L_c;  
            end       
            sum2 = sum2 + factorial(m1)/(factorial(jjj)*factorial(m1-jjj))*2^(-m1)*sum1;
        end
        sum3 = sum3 + (-1)^jj*sum2;
    end
    v_c = v_c + factorial(m1)/(factorial(j)*factorial(m1-j))*2^(-m1)*sum3;
end
for j = 1:m1    %approximation 2.2
    sum3 = zeros(N,1);
    for jj = 0:m2+j    %summation in approximation 2.2
        mu = A1/(2*t1)-1i*pi/t1-1i*(-jj+1)*pi/t1;
        sum2 = zeros(N,1);
        for jjj = 0:m1    %approximation 1.1
            sum1 = zeros(N,1);
            for jjjj = 0:m2+jjj    %summation in approximation 1.1
                theta = A2/(2*t2)-1i*pi/t2-1i*(jjjj-1)*pi/t2;
                f_c = (theta*D + mu*I-G)\ones(N,1);
                L_c = 1/(theta^2)*f_c - 1/((theta^2)*mu) + (x')/(theta*mu*(mu-r));
                sum1 = sum1 + (-1)^jjjj*L_c;

            end       
            sum2 = sum2 + factorial(m1)/(factorial(jjj)*factorial(m1-jjj))*2^(-m1)*sum1;
        end
        for jjj = 1:m1    %approximation 1.2
            sum1 = zeros(N,1);
            for jjjj = 0:m2+jjj    %summation in approximation 1.2
                theta = A2/(2*t2)-1i*pi/t2-1i*(-jjjj+1)*pi/t2;
                f_c = (theta*D + mu*I-G)\ones(N,1);
                L_c = 1/(theta^2)*f_c - 1/((theta^2)*mu) + (x')/(theta*mu*(mu-r));
                sum1 = sum1 + (-1)^jjjj*L_c;  
            end       
            sum2 = sum2 + factorial(m1)/(factorial(jjj)*factorial(m1-jjj))*2^(-m1)*sum1;
        end
        sum3 = sum3 + (-1)^jj*sum2;        
    end
    v_c = v_c + factorial(m1)/(factorial(j)*factorial(m1-j))*2^(-m1)*sum3;
end
v_c = exp((A1+A2)/2)/(4*t1*t2)*v_c;
%}

end
