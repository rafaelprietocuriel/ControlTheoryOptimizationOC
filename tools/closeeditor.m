function closeeditor()

% Close down all Editor windows.
try
    Editor=com.mathworks.mlservices.MLEditorServices;
    Editor.closeAll;
catch
    Editor=matlab.desktop.editor.getAll;
    Editor.close()
end