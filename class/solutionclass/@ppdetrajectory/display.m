function display(ppdeTrj)

if isequal(get(0,'FormatSpacing'),'compact')
    disp([inputname(1) ' =']);
    if ~isempty(ppdeTrj)
        disp('ocmatclass: ppdetrajectory')
        disp('FEM Data')
        disp(ppdeTrj.femdata)
        disp(ppdeTrj.octrajectory);
    else
        disp('empty ppdetrajectory');
    end
else
    disp(' ')
    disp([inputname(1) ' =']);
    disp(' ')
    if ~isempty(ppdeTrj)
        disp('ocmatclass: ppdetrajectory')
        disp(' ');
        disp('FEM Data')
        disp(ppdeTrj.femdata)
        disp(ppdeTrj.octrajectory);
    else
        disp('empty ppdetrajectory');
        disp(' ');
    end
end