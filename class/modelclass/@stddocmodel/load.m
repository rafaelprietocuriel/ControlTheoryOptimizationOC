function varargout=load(ocObj0,varargin)
%
% LOAD loads the data for oc model into the workspace
%
% LOAD(OCOBJ) loads the data for oc model into the workspace
%
% LOAD(OCOBJ,'WORKSPACE') loads data previously stored by the
% SAVE(OCOBJ,'WORKSPACE') command.
%
% LOAD(OCOBJ,'MODELDATA') loads model structure of the model into the
% workspace and variable name 'ocStruct'. 
%
% LOAD(OCOBJ,FORMAT) uses the format string FORMAT (see SPRINTF for
% details) for the parameter values used for the file name generation.
%
% LOAD(OCOBJ,FORMAT,FORCE) FORCE = 1 forces to overwrite the results of
% OCOBJ by the data stored from a previous session. FORCE = 0 if a data
% file exists the user is asked if s/he wants to proceed and overwrite the
% model data.
%
% LOAD(OCOBJ,FORMAT,FORCE,FN) FN provides an alternative filename.
%
% LOAD(OCOBJ,FORMAT,FORCE,IDX) filename is generated of the model
% name and the parameter values for the IDX'th parameter values. 
%
% LOAD(OCOBJ,FORMAT,FORCE,FN,SAVEDIR) SAVEDIR provides an alternative
% folder. 

format='';
force=[]; % force overwriting an already existing optimal control model
fn='';
loaddir='';

if nargin==2
    if strcmpi(varargin{1},'workspace')
        workspacename=wsname(ocObj0);
        if ~exist(workspacename,'file')
            ocmatmsg('File ''%s'' does not exist on MATLAB path.\n',workspacename)
            if nargout==1
                varargout{1}=false;
            end
            return
        end
        evalin('caller',['load(''' workspacename ''')'])
        return
    elseif strcmpi(varargin{1},'modeldata')
        if evalin('caller','exist(''ocStruct'',''var'')')==1
            while 1
                answer=input('Variable ''ocStruct'' alrerady exists. Overwrite it?  (y)/n: ','s');
                if isempty(answer)
                    % default value 'y'
                    answer='y';
                end
                if strcmpi(answer,'n')
                    return
                elseif strcmpi(answer,'y')
                    break
                end
            end
        end
        evalin('caller',['ocStruct=loadmodeldata(modelname(' inputname(1) '));'])
        return
    end
end

if isempty(ocObj0)
    ocmatmsg('Optimal control model is empty. Nothing to be loaded.\n')
    if nargout==1
        varargout{1}=false;
    end
    return
end

if nargin>=2
   format=varargin{1};
end

if nargin>=3
   force=varargin{2};
end

if nargin>=4
   fn=varargin{3};
end

if nargin>=5
   loaddir=varargin{4};
end

if isempty(force)
    force=0;
end

if isempty(fn)
    fn=filename(ocObj0,format);
elseif isnumeric(fn)
    fn=filename(ocObj0,format,fn);
elseif ischar(fn)
    tmp=parameterindex(ocObj0,fn);
    if ~isempty(tmp)
        fn=filename(ocObj0,format,tmp);
    end
end

if isempty(loaddir)
    loaddir=fullocmatfile(userdatafolder(ocObj0));
end

if isempty(strfind(fn,'.mat'))
    fn=[fn '.mat'];
end
filename_mat=fullfile(loaddir,fn);

if ~exist(filename_mat,'file')
    ocmatmsg('No data file\n''%s''\nexists.\n',filename_mat)
    if nargout==1
        varargout{1}=false;
    end
    return
end
if ~force
        flag=input(['The results in the actual model ''' inputname(1) ''' are overwritten. Proceed? y/(n) : '],'s');
        if isempty(flag)
            flag='n';
        end
        % default value 'n'
        if strcmpi(flag,'n') 
            ocmatmsg('Return without loading.\n')
            if nargout==1
                varargout{1}=false;
            end
            return
        end
end
load(filename_mat)
if ~exist('ocObj','var')
    ocmatmsg('''%s'' is not a data file.\n',filename_mat)
    if nargout==1
        varargout{1}=false;
    end
    return
end
if ~isidenticalmodel(ocObj0,ocObj)
        flag=input(['The stored model is not identical to the input model.\nReturn without loading. Proceed? y/(n) : '],'s');
        if isempty(flag)
            flag='n';
        end
        % default value 'n'
        if strcmpi(flag,'n') 
            ocmatmsg('Return without loading.\n')
            if nargout==1
                varargout{1}=false;
            end
            return
        end
end
assignin('caller',inputname(1),ocObj);
if nargout==1
    varargout{1}=true;
end
