function arcinfoStruct=arcinfo(ocObj,arcarg)

arcinfoStruct=struct([]);
if isempty(ocObj)
    return
end
if nargin==1
    arcarg=0;
end
arcinfoStruct=feval(ocObj,'ArcInfo',arcarg);