% When m=1,n=1, the code for Equation (46).
function posimn11=posimn11(rD,RD,z,omega,lambda)
posimn11=(((1./(-(rD)^2)).*sinh(sqrt(z.*g(z,omega,lambda)).*rD))+(1./rD).*sqrt(z.*g(z,omega,lambda)).*cosh(sqrt(z.*g(z,omega,lambda)).*rD)).*((1./(-(RD)^2))...
        .*cosh(sqrt(z.*g(z,omega,lambda))*RD)+(1./RD).*sqrt(z.*g(z,omega,lambda)).*sinh(sqrt(z.*g(z,omega,lambda)).*RD))-(((1./(-(rD)^2)).*cosh(sqrt(z.*g(z,omega,lambda))...
        .*rD))+(1/rD).*sqrt(z.*g(z,omega,lambda)).*sinh(sqrt(z.*g(z,omega,lambda)).*rD)).*((1/(-(RD)^2)).*sinh(sqrt(z.*g(z,omega,lambda)).*RD)+(1/RD)...
        .*sqrt(z.*g(z,omega,lambda))*cosh(sqrt(z.*g(z,omega,lambda)).*RD));
