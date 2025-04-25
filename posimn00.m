% When m=0,n=0, the code for Equation (46).
function posimn00=posimn00(rD,RD,z,alpha3,alpha4,omega,lambda)
posimn00=((1./rD).*sinh(sqrt(z.*g(z,omega,lambda)).*rD).*((1./RD).*cosh(sqrt((z-(alpha3.*RD+alpha4)).*g((z-(alpha3.*RD+alpha4)),omega,lambda)).*RD)))...
        -((1./rD).*(cosh(sqrt(z.*g(z,omega,lambda))).*rD).*((1/RD).*sinh(sqrt((z-(alpha3.*RD+alpha4)).*g((z-(alpha3.*RD+alpha4)),omega,lambda)).*RD)));