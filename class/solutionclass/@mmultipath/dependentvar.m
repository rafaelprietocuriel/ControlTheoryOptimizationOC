function out=dependentvar(ocMultiPath)

n=numberofparts(ocMultiPath);
out=cell(n,1);
for ii=1:n
    out{ii}=dependentvar(ocMultiPath.solutionclass{ii});
end