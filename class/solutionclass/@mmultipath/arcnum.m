function out=arcnum(ocMultiPath)

n=numberofparts(ocMultiPath);
out=cell(n,1);
for ii=1:n
    out{ii}=arcnum(ocMultiPath.solutionclass{ii});
end