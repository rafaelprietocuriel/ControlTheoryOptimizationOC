function display(dgObj)


if isequal(get(0,'FormatSpacing'),'compact')
    disp([inputname(1) ' =']);
    if ~isempty(dgObj)
        [parval,parvar]=parametervalue(dgObj);
        disp('ocmatclass: differentialgame')
        disp(['    modelname : ' modelname(dgObj)]);
        for ii=1:parameternum(dgObj)
            disp(['    ' parvar{ii} ' : ' num2str(parval(ii))]);
        end
    else
        disp('empty differentialgame');
    end
else
    disp(' ')
    disp([inputname(1) ' =']);
    disp(' ')
    if ~isempty(dgObj)
        [parval,parvar]=parametervalue(dgObj);
        disp('ocmatclass: differentialgame')
        disp(' ');
        disp(['    modelname : ' modelname(dgObj)]);
        disp(' ')
        for ii=1:parameternum(dgObj)
            disp(['    ' parvar{ii} ' : ' num2str(parval(ii))]);
        end
    else
        disp('empty differentialgame');
    end
    disp(' ');
end