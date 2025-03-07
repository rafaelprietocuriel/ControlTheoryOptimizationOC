function [ocStruct,filename]=loadmodeldata(modelname)
%
% LOADMODELDATA(FILENAME) load model structure from file 'FILENAME'.
%
% LOADMODELDATA(MODELNAME) load model structure from file created
% during the initialization process. This file is located at
% OCMATPATH/USERMODELFOLDER/DATA, and its name is
% modelnameModelDataStructure.mat.
%
filename='';
if ~isempty(regexp(modelname,'(.mat)\>','ONCE'))
    filename=modelname;
end
if isempty(filename)
    fullfolder=fullocmatfile(getocmatfolder('userdata','',modelname));
    filename=fullfile(fullfolder,[modelname 'ModelDataStructure.mat']);
    if ~exist(filename,'file')
        ocmatmsg('Data file ''%s'' for model ''%s'' does not exist.\n',filename,modelname)
    end
end
try
    load(filename);
catch
    ocStruct=[];
    modelfolder=fullocmatfile(getocmatfolder('specificmodel','',modelname));
    if ~exist([modelname '.ocm'],'file')
        ocmatmsg('Initialization file for ''%s'' does not exist.\nCreate initialization file and restart initialization.\n',modelname)
        return
    end
    if ~exist(modelfolder,'dir')
        ocmatmsg('Model folder does not exist.\nProcess initialization file.\n')
        return
    end
    if ~exist(modelfolder,'dir')
        ocmatmsg('Data folder does not exist.\nProcess initialization file.\n')
        return
    end
    ocmatmsg('%s\n',lasterr)
end
if ~exist('ocStruct','var')
    ocmatmsg('No variable ''ocStruct'' stored in ''%s''. Empty value returned.\nRepeat initialization process for model.\n',filename)
    ocStruct=[];
end