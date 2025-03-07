function out=dependentvar(ocgLim)

out=[state(ocgLim.ocgradtrajectory); ...
    costate(ocgLim.ocgradtrajectory);];
