function display(mComp)

if isequal(get(0,'FormatSpacing'),'compact')
    disp([inputname(1) ' =']);
    if ~isempty(mComp)
        disp('ocmatclass: modelcomposite')
        disp(mComp.mComp);
    else
        disp('empty modelcomposite');
    end
    disp('consisting of ')
else
    disp(' ')
    disp([inputname(1) ' =']);
    disp(' ')
    if ~isempty(mComp)
        disp('ocmatclass: modelcomposite')
        disp(' ');
        disp(mComp.Model);
    else
        disp('empty modelcomposite');
        disp(' ');
    end
    disp('consisting of ')
end
for ii=1:order(mComp)
    display(mComp.Model{ii})
end