function n=arcnum(dgObj)

n=[];
if isempty(dgObj)
    return
end
info=retrievemodelinformation(dgObj.Model,'arcnum');
n=info.value;