function out=modelname(ocObj)

out=[];
if isempty(ocObj)
    return
end

out=cell(1,order(ocObj));
for ii=1:order(ocObj)
    out{ii}=modelname(ocObj.Model{ii});
end
