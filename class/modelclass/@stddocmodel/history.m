function history(ocObj,openflag)
%
% HISTORY appends current history file to model specific history file
%
% HISTORY(OCOBJ) creates or opens a file '[modelname]History.m' in the
% default model data folder and appends current MATLAB history file.
%
% HISTORY(OCOBJ,'OPEN') opens the history file '[modelname]History.m' in
% the MATLAB editor.

if isempty(ocObj)
    ocmatmsg('Model is empty.\n')
    return
end

if nargin==2 && strcmpi(openflag,'open')
    edit(fullocmatfile(userdatafolder(ocObj),[modelname(ocObj) 'History.m']))
    return
end

history(modelname(ocObj),fullocmatfile(userdatafolder(ocObj)))