function out=jumpargument(ocMultiPath)

n=multiplicity(ocMultiPath);
out=cell(n,1);
for ii=1:n
    out{ii}=jumpargument(ocMultiPath.solutionclass{ii});
end