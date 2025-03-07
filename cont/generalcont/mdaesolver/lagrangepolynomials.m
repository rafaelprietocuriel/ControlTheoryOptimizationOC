function prod=lagrangepolynomials(n,rho,nr)
% calculates the Lagrangepolynomials and its integrals, rho are the
% collocation points, n is the nth Lagrangepolynomial and nr specifies if
% the integral is calculated nr=1 or not nr=0. 

prod=1;
for ii=1:length(rho)
    if (ii~=n)
        prod=conv(prod,[1 -rho(ii)])/(rho(n)-rho(ii));
    end
end
if nr==1
    prod=polyint(prod);
end
