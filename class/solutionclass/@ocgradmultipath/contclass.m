function out=contclass(ocgTrj,varargin)
out='';
if isempty(ocgTrj) || ~isfield(ocgTrj.solverinfo,'contclass')
    return
end

out=ocgTrj.solverinfo.contclass;