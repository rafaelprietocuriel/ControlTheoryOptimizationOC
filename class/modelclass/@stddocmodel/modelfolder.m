function folder=modelfolder(ocObj)
%
% 
if isempty(ocObj)
    ocmatmsg('Oc model is empty.\n')
    return
end
folder=getocmatfolder('specificmodel',modeltype(ocObj),modelname(ocObj));