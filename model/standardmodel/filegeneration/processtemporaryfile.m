function modelfilename=processtemporaryfile(temporaryfile,varStruct)
%
% PROCESSTEMPORARYFILE the temporary file was derived from a template file
% and may include file generation commands (if clauses, next loops and
% addarccases). These commands are processed and the final text is written
% into an m file of the same name like the temporary file. The temporary
% file is then deleted.

% 'Whitespace' defines the White-space characters (' \b\t') and leading
% white-space characters are not included in the proecessing of the fields.
% To keep tabulators at the beginning of a line they are removed from the
% 'Whitespace' list by setting 'Whitespace',' \b'. This allows to maintain
% the indention of the template files and inclusion of empty lines.

if verLessThan('matlab','7.6.0')
    tmpFileText=readfile2cell(temporaryfile,'%[^\n\r]','CommentStyle','//','Whitespace',' \b','BufSize',1638000);
else
    tmpFileText=readfile2cell(temporaryfile,'%[^\n\r]','CommentStyle','//','Whitespace',' \b');
end
% replace file extension tmp by m
modelfilename=regexprep(temporaryfile,'(\.tmp)\>','\.m');
try
    templatemodelfilename=regexp(regexprep(tmpFileText(1,:),'\s',''),['(?<=\=)\w*'],'match');
    modelfilename=regexprep(modelfilename,'\w*(?=(\.m))',templatemodelfilename{1});
end
% process if clauses
tmpFileText=processcommand(tmpFileText,'IFCLAUSE',[],temporaryfile);

% process addarccases
while any(strcmp(tmpFileText,'!STARTADDARCCASE!'))
    tmpFileText=processcommand(tmpFileText,'CASESWITCH',varStruct,temporaryfile);
end

% process add gneral cases
while any(strcmp(tmpFileText,'!STARTADDCASE!'))
    tmpFileText=processcommand(tmpFileText,'GENERALCASESWITCH',varStruct,temporaryfile);
end
% process for next loop
tmpFileText=processcommand(tmpFileText,'FORNEXT',varStruct,temporaryfile);

% write text into file if it does not exist
copyfile(temporaryfile,modelfilename)
% delete the temporary file
delete(temporaryfile)