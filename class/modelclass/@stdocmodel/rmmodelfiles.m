function rmmodelfiles(ocObj,varargin)

modelfiles=getmodelfilenames(ocObj.Model);

modelfilepath=fullfile(getocmatpath(),modelfolder(ocObj));
mname=modelname(ocObj);

for ii=1:length(modelfiles)
    delete(fullfile(modelfilepath,[mname modelfiles(ii).name '.m']));
end