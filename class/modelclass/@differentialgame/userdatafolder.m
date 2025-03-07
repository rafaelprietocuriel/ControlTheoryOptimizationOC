function folder=userdatafolder(dgObj)
%
% 
if isempty(dgObj)
    ocmatmsg('Oc model is empty.\n')
    return
end
folder=getocmatfolder('userdata',modeltype(dgObj),modelname(dgObj));