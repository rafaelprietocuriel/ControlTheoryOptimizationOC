function folder=userdatafolder(mmObj)
%
% 
if isempty(mmObj)
    ocmatmsg('Oc model is empty.\n')
    return
end
folder=getocmatfolder('userdata',modeltype(mmObj),modelname(mmObj));