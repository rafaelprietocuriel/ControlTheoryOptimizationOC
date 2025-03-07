function [modelfolder,datamodelfolder]=generatemodelfolder(modelname)
%
%

modelfolder=getocmatfolder('specificmodel','',modelname);
fullfolder=fullocmatfile(modelfolder);
if ~exist(fullfolder,'dir')
    answer=input([strrep(fullfolder,filesep,'\\') ' does not exist. Create it?  (y)/n: '],'s');
    if isempty(answer)
        % default value 'y'
        answer='y';
    end
    if strcmpi(answer,'n') 
        ocmatmsg('Return without generating folder %s.\n',fullfolder)
        modelfolder='';
        datamodelfolder='';
        return
    else
        mkdir(fullfolder)
    end
end
datamodelfolder=getocmatfolder('userdata','',modelname);
fullfolder=fullocmatfile(datamodelfolder);
if ~exist(fullfolder,'dir')
    answer=input([strrep(fullfolder,filesep,'\\') ' does not exist. Create it?  (y)/n: '],'s');
    if isempty(answer)
        % default value 'y'
        answer='y';
    end
    if strcmpi(answer,'n') 
        ocmatmsg('Return without generating folder %s.\n',fullfolder)
        modelfolder='';
        return
    else
        mkdir(fullfolder)
    end
end