function out=dsalvagedct(mmObj,solObj)
%
% 

out=[];
if isempty(mmObj)
    return
end
par=parametervalue(mmObj);
stageindex=stage(mmObj);
if ~ismmultipath(solObj)
    return
end

if stageindex>numberofparts(solObj)
    return
end
t=time(mmObj,solObj(stageindex),1);
depvar=dependentvar(solObj(stageindex));
ct=connectiontime(solObj);
% return optimal control value evaluated at 'depvar'
out=feval(mmObj,'DiscountedSalvagevalueDerivativeConnectionTime',t(end),depvar(:,end),par,[],ct);
