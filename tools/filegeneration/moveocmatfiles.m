function varargout=moveocmatfiles(ocStruct,modelfiles,varargin)
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
outdir='';
opt=[];
if isempty(ocStruct)
    return;
end

if nargin>=3
    opt=varargin{1};
end
if nargin>=4
    force=varargin{2};
end
if nargin>=5
    targetdir=varargin{3};
end
if nargin>=6
    outdir=varargin{4};
end
if isempty(opt)
    opt=defaultocoptions;
end
if isempty(force)
    force=0;
end
MessageDisplay=strcmp(getocoptions(opt,'INIT','MessageDisplay'),'on');

actmodelname=modelname(ocStruct);
actmodeltype=modeltype(ocStruct);
if isempty(targetdir)
    targetdir=fullocmatfile(getocmatfolder('specificmodel',actmodeltype,actmodelname));
end

% folder where the automatically generated files are placed
if isempty(outdir)
    outdir=fullocmatfile(getocmatfolder('out',actmodeltype));
end

% test if modelfiles exist in the out folder
for ii=numel(modelfiles):-1:1
    filename=modelfiles(ii).filename;
    existflag=exist(filename,'file');
    if ~existflag
        modelfiles(ii)=[];
        ocmatmsg('File %s in output folder %s does not exist.\n',filename,outdir)
    end
end

% prove if targetdir exists otherwise create
existflag=exist(targetdir,'dir');

if ~existflag
    answer=input([strrep(targetdir,'\','\\') ' does not exist. Create it?  y/(n): '],'s');
    if isempty(answer)
        % default value 'y'
        answer='y';
    end
    if strcmpi(answer,'n')
        disp('Return without moving files.')
        return
    else
        mkdir(targetdir)
    end
end

if ispc
    status=zeros(numel(modelfiles),1);
    for ii=1:numel(modelfiles)
        actualfilename=char(regexp(modelfiles(ii).filename,'\w*(?=(\.m))','match'));%[actmodelname modelfiles(ii).name '.m'];
        actualfilename=[actualfilename '.m'];
        if force
            [status(ii),sysmsg]=system(['move /Y "' fullfile(outdir,actualfilename) '" "' targetdir filesep '" ']);
        else
            answer=[];
            existflag=exist(fullfile(targetdir,actualfilename),'file');
            if existflag
                answer=input(['''' actualfilename ''' already exists in "' strrep(targetdir,'\','\\') '". Overwrite it?  y/a/(n)/q: '],'s');
                if isempty(answer)
                    % default value 'y'
                    answer='n';
                end
            end
            if strcmpi(answer,'a')
                force=1;
                answer='y';
            end
            if isempty(answer) || strcmpi(answer,'y')
                [status(ii),sysmsg]=system(['move /Y "' fullfile(outdir,actualfilename) '" "' targetdir filesep '" ']);
            elseif strcmpi(answer,'q')
                return
            end
            status(ii)=2;
            if ~strcmpi(answer,'y') && existflag
                ocmatmsg('File ''%s'' not overwritten.\n',actualfilename)
            end
        end
    end
    if status(ii)==1
        ocmatmsg('Problem with file ''%s''.\n',actualfilename)
        ocmatmsg(sysmsg)
    elseif MessageDisplay
        ocmatmsg('File ''%s'' moved to ''%s''.\n',[actmodelname modelfiles(ii).name '.m'],targetdir)
    end
end

if nargout>=1
    varargout{1}=status;
end