function save(dgObj,varargin)
%
% SAVE saves an optimal control model
%
% SAVE(OCOBJ) saves the data of an optimal control model under the standard
% data directory with the name composed of the name of the model and its
% parameter values.
%
% SAVE(OCOBJ,'WORKSPACE') stores all variables from the current MATLAB
% workspace into the model specific mat-file [modelname]Workspace, located
% in the default model data folder.
%
% SAVE(OCOBJ,FORMAT) uses the format string FORMAT (see SPRINTF for
% details) for the parameter values.
%
% SAVE(OCOBJ,FORMAT,FORCE) FORCE = 1 forces to overwrite an already
% existing file. Default value is 0.
%
% SAVE(OCOBJ,FORMAT,FORCE,FN) FN provides an alternative filename.
%
% SAVE(OCOBJ,FORMAT,FORCE,[],IDX) filename is generated of the model
% name and the parameter values for the IDX'th parameter values. 

format='';
force=[]; % force overwriting an already existing optimal control model
fn='';
savedir='';

if nargin==2 && strcmpi(varargin{1},'workspace')
    workspacename=wsname(dgObj);
    if isempty(workspacename)
        ocmatmsg('To save workspace model ''%s'' has to be non-empty.\n',inputname(1))
        return
    end
    evalin('base',['save(''' workspacename ''')']);
    return
end 

if isempty(dgObj)
    warning('Optimal control model is empty. Nothing to be stored.')
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
   if ischar(fn)
       tmp=parameterindex(dgObj,fn);
       if ~isempty(tmp)
           fn=tmp;
       end
   end
end

if nargin>=5
   savedir=varargin{4};
end

if isempty(force)
    force=0;
end

if isempty(fn)
    fn=filename(dgObj,format);
elseif isnumeric(fn)
    fn=filename(dgObj,format,fn);
end

if isempty(savedir)
    savedir=fullocmatfile(userdatafolder(dgObj));
end

filename_mat=fullfile(savedir,[fn '.mat']);
if ~exist(savedir,'dir')
    answer=input([strrep(savedir,filesep,'\\') ' does not exist. Create it?  (y)/n: '],'s');
    if isempty(answer)
        % default value 'y'
        answer='y';
    end
    if strcmpi(answer,'n') 
        ocmatmsg('Return without saving.\n')
        return
    else
        mkdir(savedir)
    end
end
if ~force
    tmp_file=fullfile(savedir,')(.tmp');
    if exist(filename_mat,'file')==2
        d=dir(filename_mat);
        save(tmp_file,'dgObj')   
        d1=dir(tmp_file);
        delete(tmp_file)
        flag=input(['Overwrite existing file from ' d.date ' of size ' num2str(d.bytes) ' with new size ' num2str(d1.bytes) '?  y/(n): '],'s');
        if isempty(flag)
            flag='n';
        end
        % default value 'n'
        if strcmpi(flag,'n') 
            ocmatmsg('Return without saving.\n')
            return
        end
    end
end

if ~ispc % adapt differences to pc and unix
    filename_mat = strrep(filename_mat, 'e+', 'e+0');
    filename_mat = strrep(filename_mat, 'e-', 'e-0');
end

save(filename_mat,'dgObj')