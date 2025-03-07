function history(dgObj,openflag)
%
% HISTORY appends current history file to model specific history file
%
% HISTORY(dgObj) creates or opens a file '[modelname]History.m' in the
% default model data folder and appends current MATLAB history file.
%
% HISTORY(dgObj,'OPEN') opens the history file '[modelname]History.m' in
% the MATLAB editor.

% if isempty(dgObj)
%     ocmatmsg('Model is empty.\n')
%     return
% end
% 
if nargin==2 && strcmpi(openflag,'open')
    edit(fullocmatfile(userdatafolder(dgObj),[modelname(dgObj) 'History.m']))
    return
end

if isempty(dgObj)
    history('General',fullocmatfile(getocmatfolder('usermodel')))
else
    history(modelname(dgObj),fullocmatfile(userdatafolder(dgObj)))
end