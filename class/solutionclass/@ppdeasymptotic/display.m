function display(ppdeAsym)

if isequal(get(0,'FormatSpacing'),'compact')
    disp([inputname(1) ' =']);
    if ~isempty(ppdeAsym)
        disp('ocmatclass: ppdeasymptotic')
        display(ppdeAsym);
    else
        disp('empty ppdeasymptotic');
    end
else
    disp(' ')
    disp([inputname(1) ' =']);
    disp(' ')
    if ~isempty(ppdeAsym)
        disp('ocmatclass: ppdeasymptotic')
        disp(' ');
        display(ppdeAsym.limitset);
        display(ppdeAsym.ppdetrajectory);
    else
        disp('empty ppdeasymptotic');
        disp(' ');
    end
end