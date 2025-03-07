function out=implicitcontrolcoordinate(dgObj,arcid)
%
out=[];
if isempty(dgObj)
    return
end
if nargin==1 || isempty(arcid)
    arcid=0;
end
if ~ischar(arcid)
    arcid=num2str(arcid);
end
info=retrievedifferentialgameinformation(dgObj.Model,'implicitnonlinearcontrolindex',arcid);
out=info.value;