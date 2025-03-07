function ocgTrj=refinetimegrid(ocgTrj,varargin)
%
% OCGTRJ=REFINETIMEGRID(OCGTRJ) refines the time grid of the
% ocgradtrajectory OCGTRJ adding new time points in the middle of every
% time interval. The state, costate and control are linearly interpolated
% at theses addtional time points.

if isempty(ocgTrj)
    return
end

t=time(ocgTrj);
tnew=sort([t t(1)+0.5*(t(1:end-1)+t(2:end))]);
ocgTrj=deval(ocgTrj,tnew);