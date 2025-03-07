function out=gradientsolver(ocgTrj)

solvr=solver(ocgTrj);

out=func2str(solvr.gradient);