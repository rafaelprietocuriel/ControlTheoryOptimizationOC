function display(docAsym)

if isequal(get(0,'FormatSpacing'),'compact')
    disp([inputname(1) ' =']);
    if ~isempty(docAsym)
        disp('ocmatclass: docasymptotic')
        disp(docAsym);
    else
        disp('empty docasymptotic');
    end
else
    disp(' ')
    disp([inputname(1) ' =']);
    disp(' ')
    if ~isempty(docAsym)
        disp('ocmatclass: docasymptotic')
        disp(' ');
        disp(docAsym.limitset);
        disp(docAsym.doctrajectory);
    else
        disp('empty docasymptotic');
        disp(' ');
    end
end