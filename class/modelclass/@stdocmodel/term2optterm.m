function optterm=term2optterm(ocObj,term,arcid,opt)
%
% TERM2OPTTERM replaces the control variables by its optimized value
%
%

if isempty(ocObj)
    optterm=[];
    return
end

if nargin==3
    opt=defaultocoptions;
end
load(ocObj,'modeldata');
optterm=term2optterm_standardmodel(ocStruct,term,arcid,opt);
