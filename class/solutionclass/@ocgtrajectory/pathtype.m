function out=pathtype(ocgTrj)

out=[];
if isfield(ocgTrj.octrajectory.solverinfo,'pathtype')
    out=ocgTrj.octrajectory.solverinfo.pathtype;
end