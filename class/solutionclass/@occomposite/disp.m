function disp(ocComp)

if isequal(get(0,'FormatSpacing'),'compact')
    if ~isempty(ocComp)
        disp('ocmatclass: occomposite')
        disp(ocComp.path);
    else
        disp('empty occomposite');
    end
else
    if ~isempty(ocComp)
        disp('ocmatclass: occomposite')
        disp(' ');
        disp(ocComp.path);
    else
        disp('empty occomposite');
        disp(' ');
    end
end