function display(ocgAsym)

if isequal(get(0,'FormatSpacing'),'compact')
    disp([inputname(1) ' =']);
    if ~isempty(ocgAsym)
        disp('ocmatclass: ocgasymptotic')
        disp(ocgAsym.limitset);
        disp(ocgAsym.ocgtrajectory);
    else
        disp('empty ocgasymptotic');
    end
else
    disp(' ')
    disp([inputname(1) ' =']);
    disp(' ')
    if ~isempty(ocgAsym)
        disp('ocmatclass: ocgasymptotic')
        disp(' ');
        disp(ocgAsym.limitset);
        disp(ocgAsym.ocgtrajectory);
    else
        disp('empty ocgasymptotic');
        disp(' ');
    end
end