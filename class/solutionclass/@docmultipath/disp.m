function disp(docMultiPath)

if isequal(get(0,'FormatSpacing'),'compact')
    if ~isempty(docMultiPath)
        disp('ocmatclass: docmultipath')
        disp(struct(docMultiPath));
    else
        disp('empty docmultipath');
    end
else
    if ~isempty(docMultiPath)
        disp('ocmatclass: docmultipath')
        disp(' ');
        disp(struct(docMultiPath));
    else
        disp('empty docmultipath');
        disp(' ');
    end
end