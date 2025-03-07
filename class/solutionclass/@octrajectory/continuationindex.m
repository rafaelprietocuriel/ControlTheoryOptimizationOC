function out=continuationindex(ocTrj)
out=[];

if isfield(ocTrj.solverinfo,'continuationindex')
    out=ocTrj.solverinfo.continuationindex;
end