function disp(ocCuv)
    
if isequal(get(0,'FormatSpacing'),'compact')
    if ~isempty(ocCuv)
        disp('ocmatclass: occurve')
        disp(struct(ocCuv));
    else
        disp('empty occurve');
    end
else
    if ~isempty(ocCuv)
        disp('ocmatclass: occurve')
        disp(' ');
        disp(struct(ocCuv));
    else
        disp('empty occurve');
    end
end