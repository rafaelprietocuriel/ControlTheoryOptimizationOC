function editinitfile(modelname,modeltype,varargin)
%
% EDITINITFILE


textfile=fullfile(getocmatpath(),getocmatfolder('initializationfile',modeltype,modelname),[modelname '.ocm']);
if ~exist(textfile,'file')
    ocmatmsg(['Initialization file for ' modelname ' does not exist.' ])
    return
end
if ispc
    try
        winopen(textfile);
    catch
        edit(textfile)
    end
else
    edit(textfile)
end