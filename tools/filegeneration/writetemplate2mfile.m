function mfilename=writetemplate2mfile(templatename,ocStruct,varTransStruct)

%fprintf('Actual file name: ''%s''\n',templatename)
mfilename=processtemporaryfile(writetemplate2temporaryfile(templatename,ocStruct,varTransStruct),varTransStruct);
