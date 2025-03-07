function [out descr]=canonicalsystemequationnum(ocObj,arcarg)

out=[];
if isempty(ocObj)
    return
end
info=retrievemodelinformation(ocObj.Model,'canonicalsystemequationnum',arcarg2arcidentifier(arcarg));
out=info.value;
descr=info.description;