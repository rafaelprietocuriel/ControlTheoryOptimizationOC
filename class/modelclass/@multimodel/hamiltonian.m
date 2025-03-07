function out=hamiltonian(mmObj,solObj,varargin)
out=[];
if isempty(mmObj)
    return
end

if ~ismmultipath(solObj)
    return
end

if numberofparts(solObj)~=numberofmodels(mmObj)
    return
end
for ii=1:numberofmodels(mmObj)
    out=[out hamiltonian(mmObj.Model{ii},solObj,1)];
end
