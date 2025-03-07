function out=arcargument(ocMultiPath)

n=numberofparts(ocMultiPath);
out=cell(n,1);
for ii=1:n
    out{ii}=arcargument(ocMultiPath.solutionclass{ii});
end