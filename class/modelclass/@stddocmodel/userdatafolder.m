function folder=userdatafolder(ocObj)
%
% 
if isempty(ocObj)
    ocmatmsg('Oc model is empty.\n')
    return
end
folder=getocmatfolder('userdata',modeltype(ocObj),modelname(ocObj));