function display(docMultiPath)

if isequal(get(0,'FormatSpacing'),'compact')
    disp([inputname(1) ' =']);
    if ~isempty(docMultiPath)
        disp('ocmatclass: docmultipath')
        disp(struct(docMultiPath));
    else
        disp('empty docmultipath');
    end
else
    disp(' ')
    disp([inputname(1) ' =']);
    disp(' ')
    if ~isempty(docMultiPath)
        disp('ocmatclass: docmultipath')
        disp(' ');
        disp(struct(docMultiPath));
    else
        disp('empty docmultipath');
        disp(' ');
    end
end