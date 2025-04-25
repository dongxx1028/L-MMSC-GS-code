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
omega=[0.005  0.235 0.465 0.695 0.825 ];
lambda=10^(-8);
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
        v(j)=v(j) + k^(M/2) * factorial(2*k+1) / (factorial(M/2-k+1) * factorial(k+1) * factorial(k) * factorial(j-k+1) * factorial(2*k-j+1));%¸Ä½ø
        
    end
    v(j)=(-1)^(M/2+j)*v(j);
end
% When omega=[0.005  0.235 0.465 0.695 0.825 ],reverse solution of Laplace space pseudo-bottomhole pressure (equation (47)) to real space ...
% pseudo-bottomhole pressure (equation (48)) 
for k=1:5
    for i=1:L
        ft(i,k)=0;
        u=log(2)/tD(i);
        for j=1:M
            z=j*u;
            % Code for numerator of equation 47)
            fz1=G1(rD,RD,z,alpha1,alpha2,alpha3,alpha4,omega(k),lambda)-S;
            % Code for denominator of equation(47)
            fm1=CD.*z.*G1(1,RD,z,alpha1,alpha2,alpha3,alpha4,omega(k),lambda)-(1+S.*CD.*z);
            % Code for equation (47)
            ft(i,k) =  ft(i,k) + v(j) *((1/z).*(fz1/fm1));
        end
        y(i,k)=ft(i,k)*u;
    end
end

% Three point method for derivative calculation (pseudo-bottomhole pressure derivative),when omega=[0.005  0.235 0.465 0.695 0.825 ];
for k=1:5
    D_y(1,k) = tD(1) * (-3 * y(1,k) + 4 * y(2,k) - y(3,k)) / (4 * tD(2) - 3 * tD(1) - tD(3));
    D_y(2,k) = tD(2) * (y(3,k) - y(1,k)) / (tD(3) - tD(1));
    for i=3:L
        D_y(i,k) = tD(i) * (y(i - 2,k) - 4 * y(i - 1,k) + 3 * y(i,k)) / (tD(i - 2) - 4 * tD(i - 1) + 3 * tD(i));
    end
end

% Double logarithmic curves of pseudo-bottomhole pressure at different omega
loglog(tD,y(:,1),'Color',[0, 0, 1], 'Linewidth',1.5)
hold on
loglog(tD,y(:,2),'Color',[0, 0.5, 0], 'Linewidth',1.5)
hold on
loglog(tD,y(:,3),'Color',[ 1, 0, 0], 'Linewidth',1.5)
hold on
loglog(tD,y(:,4),'Color',[0, 0.75, 0.75], 'Linewidth',1.5)
hold on
loglog(tD,y(:,5),'Color',[0.75, 0, 0.75], 'Linewidth',1.5)
hold on
loglog(tD,y(:,1),'Color',[0, 0, 1], 'Linewidth',1.5)
hold on
% Double logarithmic curves of pseudo-bottomhole pressure derivatives at different omega
loglog(tD,D_y(:,1),'--','Color',[0, 0, 1], 'Linewidth',1.5)
hold on
loglog(tD,D_y(:,2),'--','Color',[0, 0.5, 0], 'Linewidth',1.5)
hold on
loglog(tD,D_y(:,3),'--','Color',[ 1, 0, 0], 'Linewidth',1.5)
hold on
loglog(tD,D_y(:,4),'--','Color',[0, 0.75,0.75], 'Linewidth',1.5)
hold on
loglog(tD,D_y(:,5),'--','Color',[0.75, 0, 0.75], 'Linewidth',1.5)
hold on
loglog(tD,D_y(:,1),'--','Color',[0, 0, 1], 'Linewidth',1.5)
hold on 

axis([1, 10^12, 10^(-3), 10^4]);
hold on

hold on
xlabel('\itt_{D}','FontName','Times New Roman','FontSize',10) 
ylabel('\itP_{wD} , \itP\prime_{wD}' ,'FontName','Times New Roman','FontSize',10) 
legend('\it\omega=0.005','\it\omega=0.235','\it\omega=0.465','\it\omega=0.695','\it\omega=0.825')


