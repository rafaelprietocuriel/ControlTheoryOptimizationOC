function [coeff,p,q,D,coeff0]=cubicsolution(eq,var,ex)

% Lus2=sym('B*(lambda1*s*u^2+2*lambda1*s*u*tau-c*tau-lambda1-lambda1*tau-tau*lambda1*s)*(1-u)-c*(u+tau)*(1+tau)')
% Lus3=collect(Lus2,'u')
% [coeff,p,q,D]=cubicsolution(Lus3,'u',sym('u^3-u^2-2*u+1'));

p=[];
q=[];
if nargin==1 || isempty(var)
    var=mystr2sym(findsym(sym(eq)));
    ex=0;
end
if nargin==2
    ex=0;
end
var=mystr2sym(var);
eq=mystr2sym(eq);

coeff=repmat(sym(''),1,4);
coeff(4)=subs(eq,var,0);
fac=1;
for ii=1:3
    eq=diff(eq,var);
    fac=fac*ii;
    coeff(4-ii)=subs(eq,var,0)/fac;
end

%s=findsym(coeff(1));
%strcmp(s,var)
if coeff(1)==0
    return
end

coeff0=simplify(coeff/coeff(1));

r=coeff0(2);
s=coeff0(3);
t=coeff0(4);

p=simple(s-r^2/3);
q=simple(2*r^3/27-r*s/3+t);

D=simple((q/2)^2+(p/3)^3);

% Example
if ex
    [coeffex,pex,qex,Dex,coeff0ex]=cubicsolution('u^3-u^2+2*u+1');
    coeffex=double(coeffex);
    coeff0ex=double(coeff0ex);
    
    % Matlab solution
    a=coeff0ex(1);
    b=coeff0ex(2);
    c=coeff0ex(3);
    d=coeff0ex(4);
    msol(1)=1/6/a*(36*c*b*a-108*d*a^2-8*b^3+12*3^(1/2)*(4*c^3*a-c^2*b^2-18*c*b*a*d+27*d^2*a^2+4*d*b^3)^(1/2)*a)^(1/3)-2/3*(3*c*a-b^2)/a/(36*c*b*a-108*d*a^2-8*b^3+12*3^(1/2)*(4*c^3*a-c^2*b^2-18*c*b*a*d+27*d^2*a^2+4*d*b^3)^(1/2)*a)^(1/3)-1/3*b/a;
    msol(2)=-1/12/a*(36*c*b*a-108*d*a^2-8*b^3+12*3^(1/2)*(4*c^3*a-c^2*b^2-18*c*b*a*d+27*d^2*a^2+4*d*b^3)^(1/2)*a)^(1/3)+1/3*(3*c*a-b^2)/a/(36*c*b*a-108*d*a^2-8*b^3+12*3^(1/2)*(4*c^3*a-c^2*b^2-18*c*b*a*d+27*d^2*a^2+4*d*b^3)^(1/2)*a)^(1/3)-1/3*b/a+1/2*i*3^(1/2)*(1/6/a*(36*c*b*a-108*d*a^2-8*b^3+12*3^(1/2)*(4*c^3*a-c^2*b^2-18*c*b*a*d+27*d^2*a^2+4*d*b^3)^(1/2)*a)^(1/3)+2/3*(3*c*a-b^2)/a/(36*c*b*a-108*d*a^2-8*b^3+12*3^(1/2)*(4*c^3*a-c^2*b^2-18*c*b*a*d+27*d^2*a^2+4*d*b^3)^(1/2)*a)^(1/3));
    msol(3)=-1/12/a*(36*c*b*a-108*d*a^2-8*b^3+12*3^(1/2)*(4*c^3*a-c^2*b^2-18*c*b*a*d+27*d^2*a^2+4*d*b^3)^(1/2)*a)^(1/3)+1/3*(3*c*a-b^2)/a/(36*c*b*a-108*d*a^2-8*b^3+12*3^(1/2)*(4*c^3*a-c^2*b^2-18*c*b*a*d+27*d^2*a^2+4*d*b^3)^(1/2)*a)^(1/3)-1/3*b/a-1/2*i*3^(1/2)*(1/6/a*(36*c*b*a-108*d*a^2-8*b^3+12*3^(1/2)*(4*c^3*a-c^2*b^2-18*c*b*a*d+27*d^2*a^2+4*d*b^3)^(1/2)*a)^(1/3)+2/3*(3*c*a-b^2)/a/(36*c*b*a-108*d*a^2-8*b^3+12*3^(1/2)*(4*c^3*a-c^2*b^2-18*c*b*a*d+27*d^2*a^2+4*d*b^3)^(1/2)*a)^(1/3));
    %msol=double(solve(ex));
    Dex=double(Dex);
    qex=double(qex);
    pex=double(pex);
    if Dex>0
        u1=nthroot(-qex/2+sqrt((qex/2)^2+(pex/3)^3),3);
        v1=nthroot(-qex/2-sqrt((qex/2)^2+(pex/3)^3),3);
        sol(1)=u1+v1-coeffex(2)/3/coeffex(1);
        sol(2)=-(u1+v1)/2+((u1-v1)/2)*sqrt(-3)-coeffex(2)/3/coeffex(1);
        sol(3)=-(u1+v1)/2-((u1-v1)/2)*sqrt(-3)-coeffex(2)/3/coeffex(1);
    elseif Dex<0
        r=sqrt(-(pex/3)^3);
        phi=acos(-qex/2/r);
        %phi=2*pi-phi;
        sol(1)=2*r^(1/3)*cos(phi/3)-coeffex(2)/3/coeffex(1);
        sol(2)=2*r^(1/3)*cos(phi/3+pi*2/3)-coeffex(2)/3/coeffex(1);
        sol(3)=2*r^(1/3)*cos(phi/3+pi*4/3)-coeffex(2)/3/coeffex(1);
    else
    end
    for ii=1:length(msol)
        fprintf('Matlab solution %d  : %f (%f)\nimaginary part: %e\n',ii,msol(ii),double(subs(ex,'u',msol(ii))),imag(msol(ii)));
        fprintf('Cardano solution %d : %f (%f)\nimaginary part: %e\n',ii,sol(ii),double(subs(ex,'u',sol(ii))),imag(sol(ii)));
    end
end