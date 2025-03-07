function out=canonicalsystemdimension(ocObj,arcarg)

out=[];
if isempty(ocObj)
    return
end
if nargin==1
    arcarg=0;
end
odeinfo=retrievemodelinformation(ocObj.Model,'odedim');
aeinfo=retrievemodelinformation(ocObj.Model,'algebraicequationnum',num2str(arcarg));
out=odeinfo.value+aeinfo.value;