function out=pathtype(hocTrj)

out=[];
if isfield(hocTrj.solverinfo,'pathtype')
    out=hocTrj.solverinfo.pathtype;
end