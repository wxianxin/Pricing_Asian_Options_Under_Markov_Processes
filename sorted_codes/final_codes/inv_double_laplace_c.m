function  v_c = qianyun_inv_double_laplace_c(N,T,n,K,r,x,G)
% 6 Compute f->v
% Continuous case

m1 = 10;
m2 = 15;
I = diag(ones(1,N));
D = diag(x);
v_c = zeros(N,1);
A1 = 20;
A2 = 20;
t1 = T;    %time
t2 = K;    %strike price
l1=1;
l2=l1;

%part A
mu=A1/(2*l1*t1); 
theta=A2/(2*l2*t2); 
f_c =(theta*D + mu*I-G)\ones(N,1);
sumA=1/(theta^2)*f_c - 1/((theta^2)*mu) + (x')/(theta*mu*(mu-r));
 
%part B
sumB=zeros(N,1);
for k1=1:l2
   for e=0:m1 
   sumb=zeros(N,1);  
    for k=0:m2+e
      mu_b=A1/(2*l1*t1);
      theta_b=A2/(2*l2*t2)-1i*k1*pi/(t2*l2)-1i*k*pi/t2;
      f_c =(theta_b*D + mu_b*I-G)\ones(N,1);
      L_c = 1/(theta_b^2)*f_c - 1/((theta_b^2)*mu_b) + (x')/(theta_b*mu_b*(mu_b-r));
      sumb=sumb+2*(-1)^k*real(exp(-1i*k1*pi/l2)*L_c);
    end
   sumB=sumB+factorial(m1)/(factorial(e)*factorial(m1-e))*2^(-m1)*sumb;
   end
end

sJ=zeros(N,1);
for j1=1:l1
for e2=0:m1
sumj=zeros(N,1);  
for j=0:m2+e2    
sumC=zeros(N,1);
sumE=zeros(N,1);
for k1=1:l2
   for e =0:m1
   sumc=zeros(N,1);
   sume=zeros(N,1);
    for k=0:m2+e
      %part C
      mu_c=A1/(2*l1*t1)-1i*j1*pi/(t1*l1)-1i*j*pi/t1;
      theta_c=A2/(2*l2*t2)-1i*k1*pi/(t2*l2)-1i*k*pi/t2;
      f_c =(theta_c*D + mu_c*I-G)\ones(N,1);
      L_c = 1/(theta_c^2)*f_c - 1/((theta_c^2)*mu_c) + (x')/(theta_c*mu_c*(mu_c-r));
      sumc=sumc+(-1)^k*exp(-(1i*j1*pi/l1+1i*k1*pi/l2))*L_c; 
      %part E
      mu_e=A1/(2*l1*t1)-1i*j1*pi/(t1*l1)-1i*j*pi/t1;
      theta_e=A2/(2*l2*t2)+1i*k1*pi/(t2*l2)+1i*k*pi/t2;
      f_c =(theta_e*D + mu_e*I-G)\ones(N,1);
      L_c = 1/(theta_e^2)*f_c - 1/((theta_e^2)*mu_e) + (x')/(theta_e*mu_e*(mu_e-r));
      sume=sume+(-1)^k*exp(-(1i*j1*pi/l1-1i*k1*pi/l2))*L_c; 
    end
   sumC = sumC + factorial(m1)/(factorial(e)*factorial(m1-e))*2^(-m1)*sumc;  
   sumE = sumE + factorial(m1)/(factorial(e)*factorial(m1-e))*2^(-m1)*sume;
   end
end

%part D
mu_d=A1/(2*l1*t1)-1i*j1*pi/(t1*l1)-1i*j*pi/t1;
theta_d=A2/(2*l2*t2);
f_c =(theta_d*D + mu_d*I-G)\ones(N,1);
L_c = 1/(theta_d^2)*f_c - 1/((theta_d^2)*mu_d) + (x')/(theta_d*mu_d*(mu_d-r));
sumD=exp(-1i*j1*pi/l1)*L_c;

%Euler for j
sumj=sumj+2*(-1)^j*real(sumC+sumD+sumE);

end
sJ=sJ+factorial(m1)/(factorial(e2)*factorial(m1-e2))*2^(-m1)*sumj;
end
end

v_c=exp(A1/(2*l1)+A2/(2*l2))/(4*t1*t2*l1*l2)*(sumA+sumB+sJ);

      
      