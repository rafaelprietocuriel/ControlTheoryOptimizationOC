function out=independentvar(ocMultiPath)

n=multiplicity(ocMultiPath);
out=cell(n,1);
for ii=1:n
    out{ii}=independentvar(ocMultiPath.solutionclass{ii});
end