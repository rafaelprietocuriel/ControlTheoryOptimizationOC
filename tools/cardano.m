function [redcoeff,D]=cardano(eq,var)

redcoeff=sym([]);
D=sym([]);
if ischar(eq)
    eq=sym(eq);
end

if ischar(var)
    var=sym(var);
end

if isempty(eq)
    return
end

coeff=repmat(sym(0),1,4);

coeff(4)=subs(eq,var,0);
for ii=1:3
    eq=(eq-coeff(4-ii+1))/var;
    coeff(4-ii)=limit(eq,var,0);
end

if coeff(1)==0
    warning([inputname(1) ' is not a cubic equation in ' inputname(2) '.'])
    return
end

coeff=coeff/coeff(1);
coeff(1)=[];
redcoeff(1)=coeff(2)-coeff(1)^2/3;
redcoeff(2)=2*coeff(1)^3/27-coeff(1)*coeff(2)/3+coeff(3);
D=(redcoeff(2)/2)^2+(redcoeff(1)/3)^3;