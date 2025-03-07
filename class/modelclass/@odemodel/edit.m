function edit(ocObj)
%
% EDIT 

extension=modelextension(ocObj);
textfile=fullocmatfile(getocmatfolder('initializationfile'),[modelname(ocObj) '.' extension]);
if ~exist(textfile,'file')
    ocmatmsg('Initialization file is not in the default directory. Search in MATLAB path.\n')
    textfile=[modelname(ocObj) '.' extension];
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
