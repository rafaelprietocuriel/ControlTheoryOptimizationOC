function dynPrim=reverse(dynPrim)

if islimitcycle(dynPrim)
    dynPrim=struct(dynPrim);
    dynPrim.period=-dynPrim.period;
    dynPrim.octrajectory.x=dynPrim.octrajectory.x(end)-fliplr(dynPrim.octrajectory.x);
    dynPrim.octrajectory.y=fliplr(dynPrim.octrajectory.y);
    dynPrim.octrajectory.arcarg=fliplr(dynPrim.octrajectory.arcarg);
    dynPrim.octrajectory.arcinterval=fliplr(dynPrim.octrajectory.arcinterval)-dynPrim.octrajectory.arcinterval(end);
    dynPrim.octrajectory.arcposition=-flipud(fliplr(dynPrim.octrajectory.arcposition)-dynPrim.octrajectory.arcposition(2,end)-1);
    if ~isempty(dynPrim.octrajectory.linearization)
        dynPrim.octrajectory.linearization=inv(dynPrim.octrajectory.linearization);
    end
    dynPrim=dynprimitive(dynPrim);
else
    return
end
