% The code for Equation (44).
function value=g(z,omega,lambda)
value=omega+sqrt(3.*(1-omega).*lambda./(5.*z)).*coth(sqrt(15.*(1-omega).*z./lambda))-lambda./(5.*z);
