function savemodeldata(modelname,ocStruct)
%
%
fullfolder=fullocmatfile(getocmatfolder('userdata','',modelname));
modelfolder=fullocmatfile(getocmatfolder('specificmodel','',modelname));
filename=fullfile(fullfolder,[modelname 'ModelDataStructure.mat']);
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

save(filename,'ocStruct')