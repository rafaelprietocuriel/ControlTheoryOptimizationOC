function t=time(ocgMTrj,varargin)
% TIME returns the time argument of the ocgradtrajectory
t=[];
if isempty(ocgMTrj)
    return
end
t=ocgMTrj.argument.t;
%degr=multiplicity(ocgMTrj);

%t=reshape(t,degr.number,[]);