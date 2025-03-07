function out=canonicalmapdimension(ocObj,arcarg)

out=[];
if isempty(ocObj)
    return
end
if nargin==1
    arcarg=0;
end
odeinfo=retrievediffmodelinformation(ocObj.Model,'odedim');
aeinfo=retrievediffmodelinformation(ocObj.Model,'algebraicequationnum',num2str(arcarg));
out=odeinfo.value+aeinfo.value;