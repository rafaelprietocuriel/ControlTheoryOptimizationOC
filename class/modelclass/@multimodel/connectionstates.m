function out=connectionstates(mmObj,ocTrjMP)

out=[];
if isempty(ocTrjMP) || isempty(mmObj) || numberofmodels(mmObj)==1
    return
end

dim=statenum(mmObj);

out=zeros(dim,numberofparts(ocTrjMP)-1);
for ii=1:numberofparts(ocTrjMP)-1
    x=state(mmObj(ii),ocTrjMP(ii));
   out(:,ii)=x(:,1); 
end
