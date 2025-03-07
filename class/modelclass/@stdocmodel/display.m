function display(ocObj)


if isequal(get(0,'FormatSpacing'),'compact')
    disp([inputname(1) ' =']);
    if ~isempty(ocObj)
        [parval,parvar]=parametervalue(ocObj);
        disp('ocmatclass: stdocmodel')
        disp(['    modelname : ' modelname(ocObj)]);
        for ii=1:parameternum(ocObj)
            disp(['    ' parvar{ii} ' : ' num2str(parval(ii))]);
        end
    else
        disp('empty stdocmodel');
    end
else
    disp(' ')
    disp([inputname(1) ' =']);
    disp(' ')
    if ~isempty(ocObj)
        [parval,parvar]=parametervalue(ocObj);
        disp('ocmatclass: stdocmodel')
        disp(' ');
        disp(['    modelname : ' modelname(ocObj)]);
        disp(' ')
        for ii=1:parameternum(ocObj)
            disp(['    ' parvar{ii} ' : ' num2str(parval(ii))]);
        end
    else
        disp('empty stdocmodel');
    end
    disp(' ');
end