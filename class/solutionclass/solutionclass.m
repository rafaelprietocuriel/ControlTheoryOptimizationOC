function classname=solutionclass(solStruct)
%
% SOLUTIONCLASS determines the class name of the solution structure.
% class is not used in the sense of object oriented programming! 

classname='';

if ~isstruct(solStruct)
    ocmatmsg('Input argument ''%s'' is not a structure.',inputname(1))
    return
end
if ~isfield(solStruct,'class')
    ocmatmsg('Input argument ''%s'' has no field ''class'' to determine solutionclass.',inputname(1))
    return
end

classname=solStruct.class;