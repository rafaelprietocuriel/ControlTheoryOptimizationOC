function errest=globalerror(solObj,ocObj,varargin)

if isoctrajectory(solObj)
    solverid=solver(solObj);
    tmesh=independentvar(solObj);
    coeff=coefficient(solObj);
end