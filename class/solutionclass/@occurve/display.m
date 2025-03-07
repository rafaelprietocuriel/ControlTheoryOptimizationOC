function display(ocCuv)
    
if isequal(get(0,'FormatSpacing'),'compact')
    disp([inputname(1) ' =']);
    if ~isempty(ocCuv)
        disp('ocmatclass: occurve')
        disp(struct(ocCuv));
    else
        disp('empty occurve');
    end
else
    disp(' ')
    disp([inputname(1) ' =']);
    disp(' ')
    if ~isempty(ocCuv)
        disp('ocmatclass: occurve')
        disp(' ');
        disp(struct(ocCuv));
    else
        disp('empty occurve');
        disp(' ');
    end
end