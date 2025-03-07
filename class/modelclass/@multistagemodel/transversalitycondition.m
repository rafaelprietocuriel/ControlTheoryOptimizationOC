function out=transversalitycondition(mmObj,solObj,varargin)
%
% TRANSVERSALITYCONDITION

if isempty(mmObj)
    out=[];
    return
end
stageindex=stage(mmObj);
par=parametervalue(mmObj);
if ismmultipath(solObj)

    if stageindex>numberofparts(solObj)
        return
    end
    arcarg=arcargument(solObj(stageindex));
    arcint=arcinterval(solObj(stageindex));
    depvar=dependentvar(solObj(stageindex));
    ct=connectiontime(solObj(stageindex));
elseif isoctrajectory(solObj)
    arcarg=arcargument(solObj);
    arcint=arcinterval(solObj);
    depvar=dependentvar(solObj);
    ct=[];
end

out=feval(mmObj,'TransversalityBC',arcint(end),depvar(:,end),par,arcarg,ct);

