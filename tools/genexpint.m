function y=genexpint(x,a)
% generalized exponential integral


y=zeros(size(x));

for ii=1:length(x)
    z=x(ii);
    fun=@(t)exp(-z.*t)./(t.^a);
    y(ii)=quadgk(fun,1,inf);
end