function out=get(ocAsym,propname)

if ischar(propname)
    if isfield(ocAsym.limitset,propname) && ~isfield(ocAsym.octrajectory,propname)
        out=ocAsym.limitset.(propname);
    elseif ~isfield(ocAsym.limitset,propname) && isfield(ocAsym.octrajectory,propname)
        out=ocAsym.limitset.(octrajectory);
    elseif isfield(ocAsym.limitset,propname) && isfield(ocAsym.octrajectory,propname)
        ocmatmsg('Property name ''%s'' appears in ''limitset'' and ''octrajectory''.\n',propname)
    end
end
if isstruct(propname)
end