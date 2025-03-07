function ocTrj=equilibrium2trajectory(dynPrim,TrivialArcMeshNum,TruncationTime)
%
%
ocTrj=[];

if ~isequilibrium(dynPrim)
    return
end
ocTrj=dynPrim.octrajectory;
ocTrj.x=linspace(0,1,TrivialArcMeshNum);
ocTrj.y=repmat(ocTrj.y,1,TrivialArcMeshNum);
ocTrj.arcposition=[1;TrivialArcMeshNum];
ocTrj.timehorizon=TruncationTime;
ocTrj.arcinterval=[0 TruncationTime];
ocTrj.solverinfo=[];
