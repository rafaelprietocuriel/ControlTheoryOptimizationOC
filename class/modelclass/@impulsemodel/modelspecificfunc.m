function func=modelspecificfunc(ocObj,funcname)
%
% returns the function handle corresponding to a function build of the
% model name and the function name FUNCNAME
func=[];

if isempty(ocObj)
    return
end
if ~ischar(funcname)
    ocmaterror('Second argument is not a string.')
end
func=str2func([modelname(ocObj) funcname]);