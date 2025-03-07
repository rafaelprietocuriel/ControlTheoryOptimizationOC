function disp(ocgAsym)

if isequal(get(0,'FormatSpacing'),'compact')
    if ~isempty(ocgAsym)
        disp('ocmatclass: ocgasymptotic')
        disp(ocgAsym);
    else
        disp('empty ocgasymptotic');
    end
else
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