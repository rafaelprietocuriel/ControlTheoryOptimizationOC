function out=lagrangefunction(ocObj)

out=[];
if isempty(ocObj)
    return
end
info=retrievemodelinformation4staticoptmodel(ocObj.Model,'lagrangefunction');
out=info.value;