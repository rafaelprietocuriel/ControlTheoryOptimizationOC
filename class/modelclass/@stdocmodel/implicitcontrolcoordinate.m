function out=implicitcontrolcoordinate(ocObj,arcid)
%
out=[];
if isempty(ocObj)
    return
end
if nargin==1 || isempty(arcid)
    arcid=0;
end
if ~ischar(arcid)
    arcid=num2str(arcid);
end
info=retrievemodelinformation(ocObj.Model,'implicitnonlinearcontrolindex',arcid);
out=info.value;