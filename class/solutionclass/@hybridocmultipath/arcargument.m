function out=arcargument(ocMultiPath)

n=multiplicity(ocMultiPath);
out=cell(n,1);
for ii=1:n
    out{ii}=arcargument(ocMultiPath.solutionclass{ii});
end