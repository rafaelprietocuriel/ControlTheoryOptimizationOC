function out=entryindex(ocMultiPath)

n=multiplicity(ocMultiPath);
out=cell(n,1);
for ii=1:n
    out{ii}=entryindex(ocMultiPath.solutionclass{ii});
end