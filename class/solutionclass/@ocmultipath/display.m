function display(ocMultiPath)

if isequal(get(0,'FormatSpacing'),'compact')
    disp([inputname(1) ' =']);
    if ~isempty(ocMultiPath)
        disp('ocmatclass: ocmultipath')
        disp(struct(ocMultiPath));
    else
        disp('empty ocmultipath');
    end
else
    disp(' ')
    disp([inputname(1) ' =']);
    disp(' ')
    if ~isempty(ocMultiPath)
        disp('ocmatclass: ocmultipath')
        disp(' ');
        disp(struct(ocMultiPath));
    else
        disp('empty ocmultipath');
        disp(' ');
    end
end