function out=pathtype(docTrj)

out=[];
if isfield(docTrj.solverinfo,'pathtype')
    out=docTrj.solverinfo.pathtype;
end