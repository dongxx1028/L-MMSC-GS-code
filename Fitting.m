clear
clc
global rD RD  z alpha1 alpha2 alpha3 alpha4 omega lambda
M=8;
rD=1;
RD=1500;
alpha1=0;
alpha2=-8;
alpha3=0;
alpha4=0;
omega=0.005;
lambda=10^(-8) ;
CD=60 ;
S=2.37;
tD = zeros(1,200);
L = length(tD);


for i=1:L
    tD(i)=0.001*10^(9*i/99);
end
% Code for V(j)(equation (36))
v = zeros(1,M);
for j=1:M   
    v(j)=0;
    for k=fix((j+1)/2):1:min(j,M/2)
        v(j)=v(j) + k^(M/2) * factorial(2*k+1) / (factorial(M/2-k+1) * factorial(k+1) * factorial(k) * factorial(j-k+1) * factorial(2*k-j+1));%改进
        
    end
    v(j)=(-1)^(M/2+j)*v(j);
end

% Reverse solution of Laplace space pseudo-bottomhole pressure (equation (47)) to real space pseudo-bottomhole pressure (equation (48)) 
    for i=1:L
        ft(i)=0;
        u=log(2)/tD(i);
        for j=1:M
            z=j*u;
            % Code for numerator of equation 47)
            fz1=G1(rD,RD,z,alpha1,alpha2,alpha3,alpha4,omega,lambda)-S;
            % Code for denominator of equation(47)
            fm1=CD.*z.*G1(1,RD,z,alpha1,alpha2,alpha3,alpha4,omega,lambda)-(1+S.*CD.*z);
           % Code for equation (47)
            ft(i) =  ft(i) + v(j) *((1/z).*(fz1/fm1));
        end
        y(i)=ft(i)*u;
    end

% Double logarithmic curves of pseudo-bottomhole pressure
loglog(tD,y(:),'Color',[0, 0, 1], 'Linewidth',1.5)
hold on

axis([100, 5*10^5, 7, 10^2]);
hold on

hold on
xlabel('\itt_{D}','FontName','Times New Roman','FontSize',10) 
ylabel('\itP_{wD} , \itP\prime_{wD}' ,'FontName','Times New Roman','FontSize',10) 

filename='1.xlsx';
A=xlsread(filename);
m=A(:,1);
n=A(:,2);
plot(m,n,'r*');
hold on
