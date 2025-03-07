function out=connectiontimenumber(msObj)

out=[];
if isempty(msObj)
    return
end
info=retrievemultistagemodelinformation(msObj.stdocmodel.Model,'connectiontimenum');
out=info.value;