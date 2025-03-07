function out=controlnum(docObj)

out=[];
if isempty(docObj)
    return
end
info=retrievemodelinformation(docObj.Model,'controlnum');
out=info.value;