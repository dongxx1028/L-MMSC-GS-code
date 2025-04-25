% The code for Equation (45).
function G1=G1(rD,RD,z,alpha1,alpha2,alpha3,alpha4,omega,lambda)
G1=((exp(alpha1.*RD+alpha2).*posimn00(rD,RD,z,alpha3,alpha4,omega,lambda))+(RD.*posimn01(rD,RD,z,omega,lambda)))...
    ./((exp(alpha1.*RD+alpha2).*posimn10(1,RD,z,alpha3,alpha4,omega,lambda))+(RD.*posimn11(1,RD,z,omega,lambda)));