function disp(pdeTrj)

if isequal(get(0,'FormatSpacing'),'compact')
    if ~isempty(pdeTrj)
        disp('ocmatclass: pdetrajectory')
        disp(struct(pdeTrj));
    else
        disp('empty pdetrajectory');
    end
else
    if ~isempty(pdeTrj)
        disp('ocmatclass: pdetrajectory')
        disp(' ');
        disp(struct(pdeTrj));
    else
        disp('empty pdetrajectory');
        disp(' ');
    end
end