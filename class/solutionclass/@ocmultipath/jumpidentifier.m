function out=jumpidentifier(ocMultiPath)

n=multiplicity(ocMultiPath);
out=cell(n,1);
for ii=1:n
    out{ii}=jumpidentifier(ocMultiPath.solutionclass{ii});
end