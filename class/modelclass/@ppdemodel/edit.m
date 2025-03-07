function edit(ppdeObj)
%
% EDIT 

textfile=fullocmatfile(getocmatfolder('initializationfile'),initfilename(ppdeObj));
if ~exist(textfile,'file')
    ocmatmsg('Initialization file is not in the default directory. Search in MATLAB path.\n')
    textfile=initfilename(ppdeObj);
    if ~exist(textfile,'file')
        ocmatmsg('Initialization file does not exist on MATLAB path.\n')
        return
    end
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
