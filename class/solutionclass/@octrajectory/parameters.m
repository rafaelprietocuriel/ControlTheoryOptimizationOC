function out=parameters(ocTrj)
out=[];

if isfield(ocTrj.solverinfo,'parameters')
    out=ocTrj.solverinfo.parameters(:).';
end