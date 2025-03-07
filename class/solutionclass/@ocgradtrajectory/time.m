function t=time(ocgTrj,varargin)
% TIME returns the time argument of the ocgradtrajectory
t=[];
if isempty(ocgTrj)
    return
end
t=ocgTrj.argument.t;