function out=dependentvariable(mmObj,solObj)

out=[];
if isempty(mmObj)
    return
end

if ~ismmultipath(solObj)
    return
end
out=cell(1,numberofparts(solObj));
for ii=1:numberofparts(solObj)
    out{ii}=dependentvariable(mmObj.Model{ii},solObj(ii));
end
