function [temporaryfile,varStruct]=writetemplate2temporaryfile(templatename,ocStruct,varargin)
%
% replace all variables which are not arc dependent and write to a
% temporary file
% template files are provided by the user for the automatic file
% generation. Theses template files are pure ascii files with ending
% '.oct'.
%
varStruct=[];

if nargin>=3
    varStruct=varargin{1};
end
actmodeltype=modeltype(ocStruct);
extension=[ocStruct.basicextension 't'];
templatefile=fullocmatfile(getocmatfolder('templatefolder',actmodeltype),[actmodeltype regexprep(templatename,['(\.' extension ')\>'],'') ['.' extension]]);
temporaryfile=fullocmatfile(getocmatfolder('out'),[modelname(ocStruct) templatename '.tmp']);

% 'Whitespace' defines the White-space characters (' \b\t') and leading
% white-space characters are not included in the proecessing of the fields.
% To keep tabulators at the beginning of a line they are removed from the
% 'Whitespace' list by setting 'Whitespace',' \b'. This allows to maintain
% the indention of the template files and inclusion of empty lines.
% Defining 'CommentStyle','//' allows to include comments in template files
% that will not appear in the automatically generated files.
templateFileText=readfile2cell(templatefile,'%[^\n\r]','CommentStyle','//','Whitespace',' \b');

% replace arc independent variables in text 'templateFileText' and write it
% into a temporary file.

fid_write=fopen(temporaryfile,'w');
linenum=numel(templateFileText);
[varname,leftidx,rightidx]=findvariableintext(templateFileText,'OCMATVARIABLE');
if isempty(varStruct)
    for ii=1:linenum
        if ~isempty(leftidx{ii})
            for jj=1:numel(leftidx{ii})
                varStruct=generatevariablestruct(ocStruct,varStruct,templateFileText{ii}((leftidx{ii}(jj)+1):(rightidx{ii}(jj)-1)));
            end
        end
    end
end
% replace arc independent variables
for ii=1:linenum
    if ~isempty(varname{ii})
        newFileText=templateFileText{ii};
        for jj=1:numel(varname{ii})
            actvarname=varname{ii}{jj}(2:end-1);
            emptystr='';
            if ~varStruct.(actvarname).arcdependent
                if ~varStruct.(actvarname).multline
                    if ischar(varStruct.(actvarname).string)
                        newFileText=strrep(newFileText,[varname{ii}{jj}],varStruct.(actvarname).string);
                    elseif iscell(varStruct.(actvarname).string)
                        newFileText=strrep(newFileText,[varname{ii}{jj}],varStruct.(actvarname).string{1});
                    end
                else
                    for kk=1:numel(varStruct.(actvarname).string)
                        if kk==1
                            emptystr=strrep(newFileText,varname{ii}{jj},'');
                            newFileText=strrep(newFileText,[varname{ii}{jj}],varStruct.(actvarname).string{kk});
                        else
                            newFileText=varStruct.(actvarname).string{kk};
                        end
                        if varStruct.(actvarname).multline==1 && kk>1
                            fprintf(fid_write,'\t');
                        end
                        if kk<numel(varStruct.(actvarname).string)
                            if ~isempty(emptystr) && kk>1
                                fprintf(fid_write,repmat('\t',1,length(emptystr)));
                            end
                            fprintf(fid_write,'%s\n',newFileText);
                        end
                    end
                end
            end
        end
        if ~isempty(emptystr)
            fprintf(fid_write,repmat('\t',1,length(emptystr)));
        end
        fprintf(fid_write,'%s\n',newFileText);
    else
        fprintf(fid_write,'%s\n',templateFileText{ii});
    end
end
fclose(fid_write);