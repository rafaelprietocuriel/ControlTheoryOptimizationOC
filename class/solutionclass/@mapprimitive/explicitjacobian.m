function J=explicitjacobian(mapPrim)

if isfixpoint(mapPrim)
    J=linearization(mapPrim);
    cmapdim=size(J,2)/2;
    J=-inv(J(:,(cmapdim+1):2*cmapdim))*J(:,1:cmapdim);
else
    J=linearization(mapPrim);
    cmapdim=size(J,1);
    J=-inv(J(:,(cmapdim+1):2*cmapdim))*J(:,1:cmapdim);
    J=J^mapPrim.period;
end