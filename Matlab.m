for i=1:110;
sig(i)=20+0.1*exp(0.12*i)

i
deter1=0; 
deter2=1;
for j =1:500;
kk(j)=j*0.002;
k=kk(j);
r_0=0.1;
l_s=0;
sigma=sig(i);
count=0;
deter=0;
s=0.00000010001;
while deter==0;
count=count+1;
s=s+0.00005;
n1=0.5;
tau=100;
beta=(sigma/(1+(tau*s).^(n1)));
alpha=(s.^2*beta+k.^2).^(0.5);
A11=k*besseli(1, k*r_0);
A12=-1i*k*besseli(1, alpha*r_0);
A13=-k*besselk(1,k*r_0);
A14=-1i*k*besselk(1,alpha*r_0);
B11=(1i*k.^2*l_s*besseli(1,k*r_0))-(1i*k*besseli(0,k*r_0));
B12=(alpha.^2*l_s*besseli(1,alpha*r_0))-(alpha*besseli(0,alpha*r_0));
B13=-(1i*k.^2*l_s*besselk(1,k*r_0))-(1i*k*besselk(0,k*r_0));
B14=(alpha.^2*l_s*besselk(1,alpha*r_0))+(alpha*besselk(0,alpha*r_0));
C11=(0.5*(alpha.^2+k.^2)*besseli(0,k))-(k*besseli(1,k))-(0.5*beta*(1-k.^2)*k*besseli(1,k));
C12=(-1i*k*alpha*besseli(0,alpha))+(1i*k*besseli(1,alpha))+(0.5*1i*k*(1-k.^2)*beta*besseli(1,alpha));
C13=(0.5*(alpha.^2+k.^2)*besselk(0,k))+(k*besselk(1,k))+(0.5*beta*(1-k.^2)*k*besselk(1,k));
C14=(1i*k*alpha*besselk(0,alpha))+(1i*k*besselk(1,alpha))+(0.5*beta*1i*k*(1-k.^2)*besselk(1,alpha));
D11=2.0*1i*k.^2*besseli(1,k);
D12=(k.^2+alpha.^2)*besseli(1,alpha);
D13=-2.0*1i*k.^2*besselk(1,k);
D14=(k.^2+alpha.^2)*besselk(1,alpha);
P11=(0.5*(alpha.^2+k.^2)*besseli(0,k))-(k*besseli(1,k))-(0.5*beta*(1-k.^2)*k*besseli(1,k));
P12=(-i*k*alpha*besseli(0,alpha))+(i*k*besseli(1,alpha))+(0.5*i*k*(1-k.^2)*beta*besseli(1,alpha));
Q11=2.0*i*k.^2*besseli(1,k);
Q12=(k.^2+alpha.^2)*besseli(1,alpha);
A = [A11, A12, A13, A14
     B11, B12, B13, B14
     C11, C12, C13, C14
     D11,D12, D13, D14
];
B2=det(A);
if count>1
if B1<0 && B2>0
deter=1;
deter2=0;
deter1=1;
end
if B1>0 && B2<0
deter=1;
deter2=0;
deter1=1;
end
end
B1=B2;
if s>1
deter=1;
if deter1==0
k1(i)=k;
end
end
if s>1
deter=1;
if deter2==0 
k2(i)=k;
deter2=1;
end
end
end
if s<0.9
ss(j)=s;
end
if s>0.9
ss(j)=0;
end;
end
smax(i)=max(ss);
[M1,I1] = max(ss);
kmax(i)=kk(I1);
end