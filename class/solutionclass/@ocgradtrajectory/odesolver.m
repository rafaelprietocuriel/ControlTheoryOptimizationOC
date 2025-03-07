function out=odesolver(ocgTrj)

solvr=solver(ocgTrj);

out=func2str(solvr.ode);