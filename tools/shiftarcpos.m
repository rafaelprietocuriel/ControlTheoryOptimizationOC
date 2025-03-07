function narcpos=shiftarcpos(arcpos,arcarg,indepvar,arcint,idx)
narcpos=[];
if numel(idx)>1
    return
end
if isempty(arcpos)
    return
end
if isempty(idx) || idx==1 || idx==arcpos(2,end)
    narcpos=arcpos;
    return
end
if arcpos(2,end)<idx
    narcpos=[];
    return
end
if any(arcpos(1,:)-idx==0) || any(arcpos(2,:)-idx==0)
    narcpos=[];
    return
end
diffarcinterval=diff(arcint);
numarc=length(arcint)-1;
s=0;
for ii=1:numarc
    arcp=arcpos(1,ii):arcpos(2,ii);
    s=s(end)+(indepvar(1,arcp)-ii+1)*diffarcinterval(ii);
    t(1,arcp)=s;
end
k=find(arcpos(1,:)-idx<0,1,'last');
arcarg=[arcarg(k:end-1) arcarg(1:k) ];
tnew=[t(idx:end) t(end)+t(2:idx)];
tnew=tnew-tnew(1);
arcintn2=arcint(2:end)-t(idx);
idxn=find(arcintn2<0,1,'last');
if isempty(idx)
    idx=0;
end
arcintn2=[arcintn2(idxn+1:end-1) arcintn2(1:idxn)+t(end)];
arcintn2=[0 arcintn2 t(end)];
narcpos=find(diff(tnew)==0);
narcpos=[1 narcpos+1;narcpos length(tnew)];

for ii=1:numarc
    arcp=narcpos(1,ii):narcpos(2,ii);
    tarc=tnew(arcp);
    arcintn(ii)=tarc(end)-tarc(1);
    s(1,arcp)=ii-1+(tarc-tarc(1))/arcintn(ii);
end
arcintn=cumsum([0 arcintn]);