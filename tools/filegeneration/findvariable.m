function varStruct=findvariable(templatename,ocStruct,varStruct)
%
% REPLACEVARIABLE replace all variables in a template file that are not arc
% dependent and write to a temporary file

if nargin==2
    varStruct=[];
end
actmodeltype=modeltype(ocStruct);
extension=[ocStruct.basicextension 't'];

templatefile=fullocmatfile(getocmatfolder('templatefolder',actmodeltype),[actmodeltype regexprep(templatename,['(\.' extension ')\>'],'') ['.' extension]]);
% Defining 'CommentStyle','//' allows to include comments in template files
% that will not appear in the automatically generated files.
templateFileText=readfile2cell(templatefile,'%[^\n\r]','CommentStyle','//');
linenum=numel(templateFileText);

% variables in a template file have to consist of upper letters only and
% have to be embraced inbetween the dollar sign
[varname,leftidx,rightidx]=findvariableintext(templateFileText,'OCMATVARIABLE');
for ii=1:linenum
    if ~isempty(leftidx{ii})
        for jj=1:numel(leftidx{ii})
            % the detected variables and its values (processed in
            % generatevariablestruct) are stored in the structure
            % 'varStruct'
            varStruct=generatevariablestruct(ocStruct,varStruct,templateFileText{ii}((leftidx{ii}(jj)+1):(rightidx{ii}(jj)-1)));
        end
    end
end