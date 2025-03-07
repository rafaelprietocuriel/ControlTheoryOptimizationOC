function out=connectiontime(msObj)

out=[];
if isempty(msObj)
    return
end
info=retrievemultistagemodelinformation(msObj.stdocmodel.Model,'connectiontime');
out=info.value;