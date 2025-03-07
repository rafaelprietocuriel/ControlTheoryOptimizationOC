function out=burmeisterfischergenparameterderivative(a,b,epsilon,mu)
% OUT=BURMEISTERFISCHER(A,B) returns the value sqrt(A.^2-B.^2)-A-B, which yields
% zero iff A>=0 and B>=0 and A*B=0 (coordinate wise)

if nargin<3
    epsilon=0;
    mu=1e-4;
end
if nargin<4
    epsilon=0;
    mu=1e-4;
end
a=a(:);
b=b(:);


out=-1./sqrt(a.^2+b.^2+mu);

