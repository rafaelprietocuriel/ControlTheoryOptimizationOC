function rmmodelpath(modelname,force)
%
% RMMODELPATH(MODELNAME) removes the model and data path from the MATLAB
% path. 
%
% RMMODELPATH(MODELNAME,FORCE) If FORCE is 1 the model and data folder are
% remove from the MATLAB path, otherwise (0) the user is asked if s/he
% wants to continue.

if nargin==1
    force=0;
end
modeldatapath=fullocmatfile(getocmatfolder('userdata','',modelname));
modelpath=fullocmatfile(getocmatfolder('specificmodel','',modelname));
if force
    rmpath(modeldatapath);
    rmpath(modelpath);
    savepath
end
searchpath=path;
if ~isempty(strfind(searchpath,[modelpath pathsep]))
    while 1
        answer=input(['Remove ' strrep(modelpath,filesep,'\\') ' from MATLAB path?  y/(n): '],'s');
        if isempty(answer)
            % default value 'n'
            answer='n';
        end
        if strcmpi(answer,'n')
            break
        elseif strcmpi(answer,'y')
            rmpath(modelpath)
            savepath
            break
        end
    end
end
if ~isempty(strfind(searchpath,[modeldatapath pathsep]))
    while 1
        answer=input(['Remove ' strrep(modeldatapath,filesep,'\\') ' from MATLAB path?  y/(n): '],'s');
        if isempty(answer)
            % default value 'n'
            answer='n';
        end
        if strcmpi(answer,'n')
            break
        elseif strcmpi(answer,'y')
            rmpath(modeldatapath)
            savepath
            break
        end
    end
end
