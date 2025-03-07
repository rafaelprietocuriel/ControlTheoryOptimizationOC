function out=arcnum(ocMultiPath)

n=multiplicity(ocMultiPath);
out=zeros(1,n);
for ii=1:n
    out(ii)=arcnum(ocMultiPath.solutionclass{ii});
end