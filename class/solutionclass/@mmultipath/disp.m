function disp(ocMultiPath)

if isequal(get(0,'FormatSpacing'),'compact')
    if ~isempty(ocMultiPath)
        disp('ocmatclass: ocmultipath')
        disp(struct(ocMultiPath));
    else
        disp('empty ocmultipath');
    end
else
    if ~isempty(ocMultiPath)
        disp('ocmatclass: ocmultipath')
        disp(' ');
        disp(struct(ocMultiPath));
    else
        disp('empty ocmultipath');
        disp(' ');
    end
end