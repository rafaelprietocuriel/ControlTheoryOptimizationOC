function out=conttype(ocgTrj,varargin)
out='';
if isempty(ocgTrj) || ~isfield(ocgTrj.solverinfo,'conttype')
    return
end

out=ocgTrj.solverinfo.conttype;