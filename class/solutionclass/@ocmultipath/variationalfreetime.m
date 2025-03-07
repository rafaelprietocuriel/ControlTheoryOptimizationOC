function out=variationalfreetime(ocMultiPath)

n=multiplicity(ocMultiPath);
out=cell(n,1);
for ii=1:n
    out{ii}=variationalfreetime(ocMultiPath.solutionclass{ii});
end