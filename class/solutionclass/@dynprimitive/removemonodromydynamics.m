function dynPrim=removemonodromydynamics(dynPrim)
%
% REMOVEMONODROMYDYNAMICS removes the monodromy part in the solution
% structure
%
% DYNPRIMN=REMOVEMONODROMYDYNAMICS(DYNPRIM) if the monodromy matrix for the
% peridoic solution DYNPRIM is calculated by a BVP method and added to the
% solution structure of the octrajectory field the monodromy part is
% removed. The such cleared solution is returned in DYNPRIMN.

if ~isfield(dynPrim.octrajectory.solverinfo,'monodromycoord')
    return
end

monodromycoord=dynPrim.octrajectory.solverinfo.monodromycoord;
dynPrim.octrajectory.y(monodromycoord,:)=[];
dynPrim.octrajectory.solverinfo=rmfield(dynPrim.octrajectory.solverinfo,'monodromycoord');

switch dynPrim.octrajectory.solver
    case 'bvp4c'
        dynPrim.octrajectory.solverinfo.yp(monodromycoord,:)=[];
    case {'bvp5c','bvp6c'}
        dynPrim.octrajectory.solverinfo.yp(monodromycoord,:)=[];
        dynPrim.octrajectory.solverinfo.ypmid(monodromycoord,:)=[];
end
