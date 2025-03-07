function display(pdeTrj)

if isequal(get(0,'FormatSpacing'),'compact')
    disp([inputname(1) ' =']);
    if ~isempty(pdeTrj)
        disp('ocmatclass: pdetrajectory')
        disp('FEM Data')
        disp(pdeTrj.femdata)
        disp(pdeTrj.pdetrajectory);
    else
        disp('empty pdetrajectory');
    end
else
    disp(' ')
    disp([inputname(1) ' =']);
    disp(' ')
    if ~isempty(pdeTrj)
        disp('ocmatclass: pdetrajectory')
        disp(' ');
        disp('FEM Data')
        disp(pdeTrj.femdata)
        disp(' ');
        disp(struct(pdeTrj));
    else
        disp('empty pdetrajectory');
        disp(' ');
    end
end