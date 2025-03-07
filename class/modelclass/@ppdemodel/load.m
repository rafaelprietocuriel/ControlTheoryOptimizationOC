function varargout=load(ppdeObj0,varargin)
%
% LOAD loads the data for oc model into the workspace
%
% LOAD(PPDEOBJ) loads the data for oc model into the workspace
%
% LOAD(PPDEOBJ,'WORKSPACE') loads data previously stored by the
% SAVE(PPDEOBJ,'WORKSPACE') command.
%
% LOAD(PPDEOBJ,'MODELDATA') loads model structure of the model into the
% workspace and variable name 'ocStruct'. 
%
% LOAD(PPDEOBJ,FORMAT) uses the format string FORMAT (see SPRINTF for
% details) for the parameter values used for the file name generation.
%
% LOAD(PPDEOBJ,FORMAT,FORCE) FORCE = 1 forces to overwrite the results of
% PPDEOBJ by the data stored from a previous session. FORCE = 0 if a data
% file exists the user is asked if s/he wants to proceed and overwrite the
% model data.
%
% LOAD(PPDEOBJ,FORMAT,FORCE,FN) FN provides an alternative filename.
%
% LOAD(PPDEOBJ,FORMAT,FORCE,IDX) filename is generated of the model
% name and the parameter values for the IDX'th parameter values. 
%
% LOAD(PPDEOBJ,FORMAT,FORCE,FN,SAVEDIR) SAVEDIR provides an alternative
% folder. 

format='';
force=[]; % force overwriting an already existing optimal control model
fn='';
loaddir='';

if nargin==2
    if strcmpi(varargin{1},'workspace')
        workspacename=wsname(ppdeObj0);
        if ~exist(workspacename,'file')
            ocmatmsg('File ''%s'' does not exist on MATLAB path.\n',workspacename)
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

if isempty(ppdeObj0)
    ocmatmsg('Optimal control model is empty. Nothing to be loaded.\n')
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
    fn=filename(ppdeObj0,format);
elseif isnumeric(fn)
    fn=filename(ppdeObj0,format,fn);
end

if isempty(loaddir)
    loaddir=fullocmatfile(userdatafolder(ppdeObj0));
end

filename_mat=fullfile(loaddir,[fn '.mat']);

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

if ~isidenticalmodel(ppdeObj0,ppdeObj)
    ocmatmsg('The model stored in\n''%s''\nis not identical to the input model.\nReturn without loading.\n',filename_mat)
    return
end
assignin('caller',inputname(1),ppdeObj);
if nargout==1
    varargout{1}=true;
end
