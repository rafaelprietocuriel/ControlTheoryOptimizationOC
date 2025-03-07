function out=exogenouscoordinate(ocMultiPath)

n=numberofparts(ocMultiPath);
out=cell(n,1);
for ii=1:n
    out{ii}=objectivevaluecoordinate(ocMultiPath.solutionclass{ii});
end