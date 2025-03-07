function o=objectivevalue(ocObj,solObj,varargin)
%
%
o=[];
if isempty(ocObj) || isempty(solObj)
    return
end

sInfo=solverinfo(solObj);

if ~isempty(sInfo) && isfield(sInfo,'objectivevaluecoord') && ~isempty(sInfo.objectivevaluecoord)
    depvar=dependentvar(solObj);
    o=depvar(sInfo.objectivevaluecoord,end);
    return
end
o=integralobjectivevalue(ocObj,solObj,varargin{:});

o=o(end)+discountedsalvagevalue(ocObj,solObj,varargin{:})+impulseobjective(ocObj,solObj,varargin{:});