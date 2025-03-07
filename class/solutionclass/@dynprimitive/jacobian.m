function J=jacobian(dynPrim)

J=[];
if isequilibrium(dynPrim)
    J=linearization(dynPrim.octrajectory);
elseif isperiodic(dynPrim)
    J=linearization(dynPrim.octrajectory);
    %J=J(:,end);
    %J=reshape(J,length(J)/2,[]);
end