function display(ocComp)

if isequal(get(0,'FormatSpacing'),'compact')
    disp([inputname(1) ' =']);
    if ~isempty(ocComp)
        disp('ocmatclass: occomposite')
        disp(ocComp.path);
    else
        disp('empty occomposite');
    end
    disp('consisting of ')
else
    disp(' ')
    disp([inputname(1) ' =']);
    disp(' ')
    if ~isempty(ocComp)
        disp('ocmatclass: occomposite')
        disp(' ');
        disp(ocComp.path);
    else
        disp('empty occomposite');
        disp(' ');
    end
    disp('consisting of ')
end
for ii=1:order(ocComp)
    display(ocComp.path{ii})
end