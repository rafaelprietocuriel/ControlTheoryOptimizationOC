function addmodelpath(modelname,force)
%
% ADDMODELPATH add model and data path to MATLAB search path.
%
% ADDMODELPATH(MODELNAME) MODELNAME is a string with the name of the
% model. 
%
% ADDMODELPATH(MODELNAME,FORCE) If FORCE is 1 the model and data folder are
% added to the MATLAB path independently if the folder is already on the
% MATLAB path. If FORCE is 0 the user is asked, if the folder is not on the
% MATLAB path

if nargin==1
    force=0;
end
modeldatapath=fullocmatfile(getocmatfolder('userdata','',modelname));
modelpath=fullocmatfile(getocmatfolder('specificmodel','',modelname));

if ~exist(modeldatapath,'dir') && ~force
    if ~exist(fullocmatfile(getocmatfolder('initializationfile'),[modelname ,'.ocm']),'file')
        ocmatmsg(['Initializationfile for ''' modelname ''' does not exist.\n'])
    else
        ocmatmsg(['First initialize ''' modelname '''.\n'])
    end
    return
end
if force
    addpath(modeldatapath);
    addpath(modelpath);
    savepath
end
% add the model folder to the MATLAB path
searchpath=path;
if isempty(strfind(searchpath,[modelpath pathsep]))
    while 1
        answer=input([strrep(modelpath,filesep,'\\') ' is not on MATLAB path. Add it?  (y)/n: '],'s');
        if isempty(answer)
            % default value 'y'
            answer='y';
        end
        if strcmpi(answer,'n')
            break
        elseif strcmpi(answer,'y')
            addpath(modelpath)
            savepath
            break
        end
    end
end
if isempty(strfind(searchpath,[modeldatapath pathsep]))
    while 1
        answer=input([strrep(modeldatapath,filesep,'\\') ' is not on MATLAB path. Add it?  (y)/n: '],'s');
        if isempty(answer)
            % default value 'y'
            answer='y';
        end
        if strcmpi(answer,'n')
            break
        elseif strcmpi(answer,'y')
            addpath(modeldatapath)
            savepath
            break
        end
    end
end
