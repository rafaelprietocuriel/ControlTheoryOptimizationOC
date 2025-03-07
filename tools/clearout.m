function clearout(ocStruct,force)
%
% CLEAROUT removes model files from the standard output folder
%
% CLEAROUT removes all files from the standard output folder.
%
% CLEAROUT(OCSTRUCT) removes the model files of OCSTRUCT from the
% standard output folder. If OCSTRUCT is empty all files are removed.
%
% CLEAROUT(OCSTRUCT,FORCE) FORCE=1 removes files without further query.

if nargin==0
    ocStruct=[];
    force=0;
end
if nargin==1
    force=0;
end
if ischar(ocStruct)
    ocStruct=loadmodeldata(ocStruct);
end
actmodelname=modelname(ocStruct);

% folder where the automatically generated files are placed
outdir=fullocmatfile(getocmatfolder('out'));

mfilesinout=dir(fullfile(outdir,[actmodelname '*.m']));
tmpfilesinout=dir(fullfile(outdir,[actmodelname '*.tmp']));
zipfilesinout=dir(fullfile(outdir,[actmodelname '*.zip']));
if isempty(mfilesinout) && isempty(tmpfilesinout) && isempty(zipfilesinout)
    fprintf('No files in the ''OCMat'' output folder.\n')
    return
end
if force
    fclose('all');
    delete(fullfile(outdir,[actmodelname '*.m']));
else
    fids=fopen('all');
    for ii=1:length(fids)
        filename=fopen(fids(ii));
        answer=input(['\nDo you want to close the file ''' strrep(filename,'\','\\') '''?  (y)/n: '],'s');
        if isempty(answer)
            % default value 'y'
            answer='y';
        end
        if strcmpi(answer,'y')
            status=fclose(fids(ii));
            if ~status
                fprintf([strrep(filename,'\','\\') ' was closed successfully.\n'])
            else
                fprintf([strrep(filename,'\','\\') ' could not be closed.\n'])
            end
        end
    end
    answer=input(['Do you want to delete the files ' [actmodelname '*.m/tmp/zip'] ' from the ''OCMat'' output folder?  y/(n): '],'s');
    if isempty(answer)
        % default value 'n
        answer='n';
    end
    if strcmpi(answer,'n')
        fprintf('\nNo files deleted.')
        return
    else
        delete(fullfile(outdir,[actmodelname '*.m']));
        delete(fullfile(outdir,[actmodelname '*.tmp']));
        delete(fullfile(outdir,[actmodelname '*.zip']));
    end
end
