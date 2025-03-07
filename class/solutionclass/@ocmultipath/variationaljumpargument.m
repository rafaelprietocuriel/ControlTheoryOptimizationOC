function out=variationaljumpargument(ocMultiPath)

n=multiplicity(ocMultiPath);
out=cell(n,1);
for ii=1:n
    out{ii}=variationaljumpargument(ocMultiPath.solutionclass{ii});
end