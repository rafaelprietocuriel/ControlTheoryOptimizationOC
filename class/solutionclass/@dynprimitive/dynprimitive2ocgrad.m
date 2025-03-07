function ocgLim=dynprimitive2ocgrad(dynPrim)

ocgLim=struct(dynPrim);
ocgLim.linearization=linearization(ocgLim.octrajectory);
ocgLim.ocgradtrajectory=octrajectory2ocgrad(ocgLim.octrajectory);
ocgLim=rmfield(ocgLim,'octrajectory');
ocgLim=ocgradlimitset(ocgLim);