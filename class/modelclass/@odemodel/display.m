function display(odeObj)


if isequal(get(0,'FormatSpacing'),'compact')
    disp([inputname(1) ' =']);
    if ~isempty(odeObj)
        [parval,parvar]=parametervalue(odeObj);
        disp('ocmatclass: odemodel')
        disp(['    modelname : ' modelname(odeObj)]);
        for ii=1:parameternum(odeObj)
            disp(['    ' parvar{ii} ' : ' num2str(parval(ii))]);
        end
    else
        disp('empty odemodel');
    end
else
    disp(' ')
    disp([inputname(1) ' =']);
    disp(' ')
    if ~isempty(odeObj)
        [parval,parvar]=parametervalue(odeObj);
        disp('ocmatclass: odemodel')
        disp(' ');
        disp(['    modelname : ' modelname(odeObj)]);
        disp(' ')
        for ii=1:parameternum(odeObj)
            disp(['    ' parvar{ii} ' : ' num2str(parval(ii))]);
        end
    else
        disp('empty odemodel');
    end
    disp(' ');
end