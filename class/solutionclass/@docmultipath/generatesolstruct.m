function sol=generatesolstruct(docMultiPath,solvername,varargin)

nummult=multiplicity(docMultiPath);
sol.x=independentvar(docMultiPath.solutionclass{1});
sol.y=dependentvar(docMultiPath.solutionclass{1});
sol.arcarg=arcargument(docMultiPath.solutionclass{1});
sol.arcinterval=arcinterval(docMultiPath.solutionclass{1});
sol.parameters=sol.arcinterval(2:arcnum(docMultiPath.solutionclass{1}));
x0=initialtime(docMultiPath);
for ii=2:nummult
    sol.x=[sol.x independentvar(docMultiPath.solutionclass{ii})-x0+sol.x(end)];
    sol.y=[sol.y dependentvar(docMultiPath.solutionclass{ii})];
    sol.arcarg=[sol.arcarg arcargument(docMultiPath.solutionclass{ii})];
    actarcinterval=arcinterval(docMultiPath.solutionclass{ii});
    sol.arcinterval=[sol.arcinterval actarcinterval];
    sol.parameters=[sol.parameters actarcinterval(2:arcnum(docMultiPath.solutionclass{ii}))];
end
sol.x0=x0;
sol.solver=solvername;

sol.solverinfo.tangent=[];
sol.solverinfo.coeff=[];
