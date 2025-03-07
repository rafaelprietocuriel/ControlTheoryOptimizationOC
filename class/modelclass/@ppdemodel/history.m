function history(ocObj)
%
% appends actual history file to model specific history file

if isempty(ocObj)
    ocmatmsg('Oc model is empty.\n')
    return
end
modelhistoryfile=fullocmatfile(userdatafolder(ocObj),[modelname(ocObj) 'History.m']);
systemhistoryfile=fullfile(prefdir,'history.m');% prefdir ... Directory containing preferences, history, and layout files

fidread=fopen(systemhistoryfile,'r+'); 
fidwrite = fopen(modelhistoryfile,'a+');
if fidread==-1
    ocmatmsg('Cannot open history file.')
    return
end
if fidwrite==-1
    ocmatmsg('Cannot open file to write %s.',modelhistoryfile)
    return
end
while 1
    tline = fgetl(fidread);
    if ~ischar(tline)
        break
    end
    fprintf(fidwrite, '%s\n',tline);
end
fclose(fidread);
fclose(fidwrite);