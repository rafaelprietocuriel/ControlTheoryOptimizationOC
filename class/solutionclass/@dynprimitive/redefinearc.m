function dynPrim=redefinearc(dynPrim,newposition,arcid)

if isperiodic(dynPrim)
    dynPrim.octrajectory=redefinearc(dynPrim.octrajectory,newposition,arcid);
end