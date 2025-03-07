function out=connectiontime(mmObj,ocTrjMP)

out=[];
if isempty(ocTrjMP) || isempty(mmObj) || numberofmodels(mmObj)==1
    return
end

out=zeros(1,numberofparts(ocTrjMP)-1);
for ii=1:numberofparts(ocTrjMP)-1
    arcint=arcinterval(ocTrjMP(ii));
    out(ii)=arcint(end);
end
