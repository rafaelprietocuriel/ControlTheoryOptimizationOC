function tmesh1_2=makefinemesh(tmesh)
global OCMATCONT
% returns the refined mesh (add new mesh point into the mid of the internal interval) based on tmesh
idx=find(diff(tmesh)==0);
arcposition=[1 idx+1;idx length(tmesh)];
N=diff(arcposition)+1;
tmesh1_2=zeros(1,2*length(tmesh)-1);
for ii=1:OCMATCONT.arcnumber
    t1=tmesh(arcposition(1,ii):arcposition(2,ii));
    tmesh1_2(1:2:2*N(ii))=t1;
    tmesh1_2(2:2:2*N(ii)-1)=t1(1:N(ii)-1)+diff(t1)/2;
end
