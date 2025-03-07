function ocTrj=composite(ocTrj1,ocTrj2)

t1=absolutetime(ocTrj1);
t2=absolutetime(ocTrj2);
if t1(end)==t2(1)
    if ocTrj1.arcarg(end)==ocTrj2.arcarg(1)
        ocTrj.arcarg=[ocTrj1.arcarg(1:end-1) ocTrj2.arcarg];
        ocTrj.x=[ocTrj1.x(1:end-1) ocTrj2.x];
        t=[t1(1:end-1) t2];
        ocTrj.y=[ocTrj1.y(:,1:end-1) ocTrj2.y];
        ocTrj.timehorizon=ocTrj1.timehorizon+ocTrj2.timehorizon;
    else
        ocTrj.arcarg=[ocTrj1.arcarg ocTrj2.arcarg];
        ocTrj.x=[ocTrj1.x ocTrj2.x];
        t=[t1 t2];
        ocTrj.y=[ocTrj1.y ocTrj2.y];
        ocTrj.timehorizon=ocTrj1.timehorizon+ocTrj2.timehorizon;
    end
else
    ocTrj=octrajectory([]);
    return
end
idx=find(diff(ocTrj.x)==0);
ocTrj.arcposition=[1 idx+1;idx length(ocTrj.x)];

[ocTrj.x,ocTrj.arcinterval]=normalizetime(t,ocTrj.arcposition);


ocTrj=octrajectory(ocTrj);

function [s arcinterval]=normalizetime(t,arcposition)

numarc=size(arcposition,2);
arcinterval=zeros(1,numarc+1);
s=t;
for ii=1:numarc
    arcidx=arcposition(ii,1):arcposition(ii,2);
    tidx=t(arcidx);
    arcinterval(ii+1)=tidx(end)-tidx(1);
    s(arcidx)=(ii-1)+(tidx-tidx(1))/arcinterval(ii+1);
end
arcinterval=cumsum(arcinterval);