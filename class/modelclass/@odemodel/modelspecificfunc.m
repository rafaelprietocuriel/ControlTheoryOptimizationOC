function func=modelspecificfunc(odeObj,funcname)
%
% returns the function handle corresponding to a function build of the
% model name and the function name FUNCNAME
func=[];

if isempty(odeObj)
    return
end
if ~ischar(funcname)
    ocmaterror('Second argument is not a string.')
end
func=str2func([modelname(odeObj) funcname]);