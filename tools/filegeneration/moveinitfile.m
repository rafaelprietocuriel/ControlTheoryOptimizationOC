function varargout=moveinitfile(modelname,varargin)
%
% MOVEOCMATFILES moves the generated files to the model folder
%
% MOVEOCMATFILES(OCOBJ,MODELFILES) moves the files 'MODELFILES' from the
% standard output folder to the standard model folder.
%
% STATUS = MOVEOCMATFILES(OCOBJ,MODELFILES,OPT,FORCE,DIR) The resulting
% status and standard output from the operating system are returned.


force=[]; % force overwriting files
targetdir='';
actmodeltype='';
if isempty(modelname)
    return;
end

if nargin>=2
    actmodeltype=varargin{1};
end
if nargin>=3
    force=varargin{2};
end
if nargin>=4
    targetdir=varargin{3};
end
if isempty(actmodeltype)
    actmodeltype='ocm';
end
if isempty(force)
    force=0;
end

initializationfile=[modelname '.' actmodeltype];
if isempty(targetdir)
    targetdir=fullocmatfile(getocmatfolder('initializationfile'));
end

localfile=fullfile('./',initializationfile);
if ~exist(localfile,'file')
    ocmatmsg('''%s'' does not exist in actual folder.\n',initializationfile)
end
% test if modelfiles exist in the out folder
existflag=exist(fullfile(targetdir,initializationfile),'file');
if existflag && ~force
    answer=input(['''' initializationfile ''' already exists in "' strrep(targetdir,'\','\\') '". Overwrite it?  y/(n): '],'s');
    if isempty(answer)
        % default value 'y'
        answer='n';
    end
    if strcmpi(answer,'a')
        answer='y';
    end
    if isempty(answer) || strcmpi(answer,'y')
        [status,sysmsg]=system(['move /Y "' localfile '" "' targetdir filesep '" ']);
    else
        return
    end
    %status=2;
    if ~strcmpi(answer,'y') && existflag
        ocmatmsg('File ''%s'' not overwritten.\n',initializationfile)
    end
else
   [status,sysmsg]=system(['move /Y "' localfile '" "' targetdir filesep '" ']); 
end

if status==1
    ocmatmsg('Problem with file ''%s''.\n',initializationfile)
    ocmatmsg(sysmsg)
else
    ocmatmsg('File ''%s'' moved to ''%s''.\n',initializationfile,targetdir)
end
if nargout>=1
    varargout{1}=status;
end