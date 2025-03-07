function [out,zeroidx,nonzeroidx]=burmeisterfischergenderivative(a,b,epsilon,mu)
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
zeroidx=find(abs(a)+abs(b)<=epsilon);
%nonzeroidx=setdiff(1:length(a),zeroidx);

nonzeroidx=1:length(a);
if ~isempty(zeroidx)
    nonzeroidx(zeroidx)=[];
end

out(nonzeroidx,1:2)=[a(nonzeroidx)./sqrt(a(nonzeroidx).^2+b(nonzeroidx).^2+mu)-1,b(nonzeroidx)./sqrt(a(nonzeroidx).^2+b(nonzeroidx).^2+mu)-1];
if ~isempty(zeroidx)
    out(zeroidx,1:2)=myrepmat([0 -1],length(zeroidx),1,1,2);
end
