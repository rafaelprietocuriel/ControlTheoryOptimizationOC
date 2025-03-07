function out=hamiltonian(ocgTrj,varargin)
out=[];
if isempty(ocgTrj)
    return
end

func=str2func([modelname(ocgTrj) 'LocalHamiltonian']);

out=func(time(ocgTrj),state(ocgTrj),control(ocgTrj),costate(ocgTrj),modelparameter(ocgTrj));