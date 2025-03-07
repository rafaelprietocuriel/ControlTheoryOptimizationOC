function B=myrepmat(A,M,N,m,n)
siz = [M N];
mind = (1:m)';
nind = (1:n)';
mind = mind(:,ones(1,siz(1)));
nind = nind(:,ones(1,siz(2)));
B = A(mind,nind);
