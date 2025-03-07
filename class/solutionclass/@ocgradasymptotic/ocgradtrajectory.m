function out=ocgradtrajectory(ocgAsym)

if isempty(ocgAsym)
    out=ocgradtrajectory();
    return
end

out=ocgAsym.ocgradtrajectory;