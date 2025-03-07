function out=continuationparameter(ocTrj)
out=[];

if isfield(ocTrj.solverinfo,'continuationparameter')
    out=ocTrj.solverinfo.continuationparameter;
end