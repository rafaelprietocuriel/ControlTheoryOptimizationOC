function varargout=rmmodel(modelname,force)
%
% RMMODEL removes a model from OCMat
%
% RMMODEL(MODELNAME) deletes the model files of model MODELNAME and removes
% the model folders. Additionally the folders are removed from the MATLAB
% path
%
% RMMODEL(MODELNAME,FORCE) FORCE=1: forces to remove the folders
% without asking.
% FORCE=0: the user is asked if s/he wants to continue.
%
% [STAT,MESS,ID]=RMMODEL(...) Here, STAT is 1 for success and is 0 for
% error, and message MESS, messageid ID.

if nargin==1
    force=0;
end
if ~isnumeric(force) || numel(force)>1
    ocmaterror('Invalid second data argument')
end

if ~force
    while 1
        accept=lower(input('All model files will be removed. Continue: y/(n) : ','s'));
        if isempty(accept)
            accept='n';
        end
        if strcmpi(accept,'y')
            break
        elseif strcmp(accept,'n')
            ocmatmsg(['Model ''' modelname ''' has NOT been removed.\n'])
            return
        end
    end
end
rmmodelpath(modelname);
modelfolder=fullocmatfile(getocmatfolder('specificmodel','',modelname));
if exist(modelfolder,'dir')
    [stat,mess,id]=rmdir(modelfolder,'s');
else
    ocmatmsg(['Model directory for ''' modelname ''' does not exist.\n'])
    return
end

if stat
    ocmatmsg(['Model ''' modelname ''' has sucessfully been removed.\n'])
end
if nargout>=1
    varargout{1}=stat;
end
if nargout>=2
    varargout{2}=mess;
end
if nargout>=3
    varargout{3}=id;
end