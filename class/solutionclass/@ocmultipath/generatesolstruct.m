function sol=generatesolstruct(ocMultiPath)

nummult=multiplicity(ocMultiPath);
sol.x=independentvar(ocMultiPath.solutionclass{1});
sol.y=dependentvar(ocMultiPath.solutionclass{1});
sol.arcarg=arcargument(ocMultiPath.solutionclass{1});
sol.arcinterval=arcinterval(ocMultiPath.solutionclass{1});
sol.parameters=parameters(ocMultiPath.solutionclass{1});
numcontpar=length(continuationparameter(ocMultiPath.solutionclass{1}));
sol.parameters(end-numcontpar+1:end)=[];
x0=initialtime(ocMultiPath);
for ii=2:nummult
    sol.x=[sol.x independentvar(ocMultiPath.solutionclass{ii})-x0+sol.x(end)];
    sol.y=[sol.y dependentvar(ocMultiPath.solutionclass{ii})];
    sol.arcarg=[sol.arcarg arcargument(ocMultiPath.solutionclass{ii})];
    actarcinterval=arcinterval(ocMultiPath.solutionclass{ii});
    sol.arcinterval=[sol.arcinterval actarcinterval];
    freepar=parameters(ocMultiPath.solutionclass{ii});
    numcontpar=length(continuationparameter(ocMultiPath.solutionclass{ii}));
    freepar(end-numcontpar+1:end)=[];
    sol.parameters=[sol.parameters freepar];
end
sol.x0=x0;
sol.solver='';

sol.solverinfo.tangent=[];
sol.solverinfo.coeff=[];
sol.solverinfo.tmesh=[];
