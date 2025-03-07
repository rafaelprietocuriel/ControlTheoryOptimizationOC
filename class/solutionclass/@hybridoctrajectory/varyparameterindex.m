function out=varyparameterindex(ocTrj)

out=[];
if isfield(ocTrj.solverinfo,'varyparameterindex')
    out=ocTrj.solverinfo.varyparameterindex;
end