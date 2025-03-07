function folder=userdatafolder(odeObj)
%
% 
if isempty(odeObj)
    ocmatmsg('Oc model is empty.\n')
    return
end
folder=getocmatfolder('userdata',modeltype(odeObj),modelname(odeObj));