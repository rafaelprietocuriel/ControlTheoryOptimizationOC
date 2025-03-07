function display(ppdeObj)


if isequal(get(0,'FormatSpacing'),'compact')
    disp([inputname(1) ' =']);
    if ~isempty(ppdeObj)
        [parval,parvar]=parametervalue(ppdeObj);
        disp('ocmatclass: ppdemodel')
        disp(['    modelname : ' modelname(ppdeObj)]);
        for ii=1:parameternum(ppdeObj)
            disp(['    ' parvar{ii} ' : ' num2str(parval(ii))]);
        end
    else
        disp('empty ppdemodel');
    end
else
    disp(' ')
    disp([inputname(1) ' =']);
    disp(' ')
    if ~isempty(ppdeObj)
        [parval,parvar]=parametervalue(ppdeObj);
        disp('ocmatclass: ppdemodel')
        disp(' ');
        disp(['    modelname : ' modelname(ppdeObj)]);
        disp(' ')
        for ii=1:parameternum(ppdeObj)
            disp(['    ' parvar{ii} ' : ' num2str(parval(ii))]);
        end
    else
        disp('empty ppdemodel');
    end
    disp(' ');
end