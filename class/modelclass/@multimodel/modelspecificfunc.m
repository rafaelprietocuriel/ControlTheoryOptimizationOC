function func=modelspecificfunc(mmObj,funcname)
%
% returns the function handle corresponding to a function build of the
% model name and the function name FUNCNAME
func=[];

if isempty(mmObj)
    return
end
if ~ischar(funcname)
    ocmaterror('Second argument is not a string.')
end
mname=submodelname(mmObj);
nummod=numberofmodels(mmObj);
for ii=1:nummod
    func{ii}=str2func([mname{ii} funcname]);
end