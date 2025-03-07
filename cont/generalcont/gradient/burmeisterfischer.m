function out=burmeisterfischer(a,b,mu)
% OUT=BURMEISTERFISCHER(A,B) returns the value sqrt(A.^2-B.^2)-A-B, which yields
% zero iff A>=0 and B>=0 and A*B=0 (coordinate wise)

if nargin==2
    mu=1e-4;
end


out=sqrt(a.^2+b.^2+mu)-a-b;