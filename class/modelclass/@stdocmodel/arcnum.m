function n=arcnum(ocObj)

n=[];
if isempty(ocObj)
    return
end
info=retrievemodelinformation(ocObj.Model,'arcnum');
n=info.value;