function out=discountrate(dgObj)
%
% DISCOUNTRATE value of the discountrate

if isempty(dgObj)
    return
end
n=playernum(dgObj);
out=zeros(n,1);

for ii=1:n
    out(ii)=subsparametervalue(dgObj,dgObj.Model.objective.(['player' num2str(ii)]).integral.discountrate);
end
%info=retrievemodelinformation(dgObj.Model,'discountrate');
%var=info.value;