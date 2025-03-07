function [out descr]=canonicalmapequationnum(ocObj,arcarg)

out=[];
if isempty(ocObj)
    return
end
info=retrievediffmodelinformation(ocObj.Model,'canonicalsystemequationnum',arcarg2arcidentifier(arcarg));
out=info.value;
descr=info.description;