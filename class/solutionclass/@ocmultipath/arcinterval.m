function out=arcinterval(ocMultiPath)

n=multiplicity(ocMultiPath);
out=cell(n,1);
for ii=1:n
    out{ii}=arcinterval(ocMultiPath.solutionclass{ii});
end