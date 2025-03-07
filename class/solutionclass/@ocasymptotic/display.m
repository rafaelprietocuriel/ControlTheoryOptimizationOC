function display(ocAsym)

if isequal(get(0,'FormatSpacing'),'compact')
    disp([inputname(1) ' =']);
    if ~isempty(ocAsym)
        disp('ocmatclass: ocasymptotic')
        disp(ocAsym);
    else
        disp('empty ocasymptotic');
    end
else
    disp(' ')
    disp([inputname(1) ' =']);
    disp(' ')
    if ~isempty(ocAsym)
        disp('ocmatclass: ocasymptotic')
        disp(' ');
        disp(ocAsym.limitset);
        disp(ocAsym.octrajectory);
    else
        disp('empty ocasymptotic');
        disp(' ');
    end
end