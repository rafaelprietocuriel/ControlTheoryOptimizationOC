function idx=rel2absidx(relsize,pos0,abssize)
%
% returns the absolute indices of a matrix A with size 'relsize', where the
% upper left point is located at pos0 within a matrix B of size 'abssize'

% example:
% A=rand(20,30);B=zeros(3,4);idx=rel2absidx(size(B),[3
% 2],size(A));B(:)=A(idx)

tmp=(0:relsize(2)-1)*abssize(1);
tmp=tmp(ones(relsize(1),1),:);
tmp2=(1:relsize(1)).';
tmp2=tmp2(:,ones(1,relsize(2)));
idx=(pos0(2)-1)*abssize(1)+pos0(1)-1+tmp+tmp2;
idx=idx(:).';