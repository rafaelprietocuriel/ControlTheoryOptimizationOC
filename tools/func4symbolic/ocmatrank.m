function idx=ocmatrank(F,v1,v2)

idx=[];

if verLessThan('symbolic','8')
    F=sym(F);
    v1=sym(v1);
    v2=sym(v2);
else
    F=str2sym(F);
    v1=str2sym(v1);
    v2=str2sym(v2);
end
if isempty(F) || isempty(v1) || isempty(v2)
    return
end
numvar=length(v2);

H=jacobian(jacobian(F,v1),v2);
c=nchoosek(1:length(v1),numvar);
for ii=1:size(c,1)
    Hpart=H(c(ii,:),:);
    if rank(Hpart)==numvar
        idx=c(ii,:);
        return
    end
end