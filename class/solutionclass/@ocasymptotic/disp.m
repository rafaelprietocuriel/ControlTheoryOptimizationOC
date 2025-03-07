function disp(ocAsym)

if isequal(get(0,'FormatSpacing'),'compact')
    if ~isempty(ocAsym)
        disp('ocmatclass: ocasymptotic')
        disp(ocAsym.limitset);
        disp(ocAsym.octrajectory);
    else
        disp('empty ocasymptotic');
    end
else
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