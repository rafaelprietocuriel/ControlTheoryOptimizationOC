function ocLSR=reversetime(ocLS)
%
% 

ocLSR=ocLS;
if isequilibrium(ocLS)
    return
end

ocLSR.octrajectory=reversetime(ocLS.octrajectory);
ocLSR.period=-ocLSR.period;