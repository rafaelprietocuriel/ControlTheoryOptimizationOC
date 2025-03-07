function out=connectiontime(ocMultiPath)

out=[];
if isempty(ocMultiPath) || numberofparts(ocMultiPath)==1
    return
end

arcint=arcinterval(ocMultiPath);

out=zeros(1,numberofparts(ocMultiPath)-1);
for ii=1:numberofparts(ocMultiPath)-1
    out(ii)=arcint{ii}(end);
end