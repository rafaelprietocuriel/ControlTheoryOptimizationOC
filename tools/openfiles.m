function OpenF=openfiles()
% return a structure consisting of the actual working directory and the
% open editor files
try
    Editor=com.mathworks.mlservices.MLEditorServices;
    openfilenames=cellstr(char(Editor.builtinGetOpenDocumentNames));
catch
    Editor=matlab.desktop.editor.getAll;
    openfilenames=cell(length(Editor),1);
    for ii=1:length(Editor)
        openfilenames{ii}=Editor(ii).Filename;
    end
end
OpenF(length(openfilenames)).localdirectory=pwd;
for ii=1:length(openfilenames)
    [pathstr,name,ext]=fileparts(openfilenames{ii});
    OpenF(ii).path=pathstr;
    OpenF(ii).filename=[name ext];
    OpenF(ii).localdirectory=OpenF(length(openfilenames)).localdirectory;
end
