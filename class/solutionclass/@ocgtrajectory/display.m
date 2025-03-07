function display(ocgTrj)

if isequal(get(0,'FormatSpacing'),'compact')
    disp([inputname(1) ' =']);
    if ~isempty(ocgTrj)
        disp('ocmatclass: ocgtrajectory')
        disp(struct(ocgTrj.octrajectory));
        disp('OdeNum:')
        disp(ocgTrj.odenum);
    else
        disp('empty ocgtrajectory');
    end
else
    disp(' ')
    disp([inputname(1) ' =']);
    disp(' ')
    if ~isempty(ocgTrj)
        disp('ocmatclass: ocgtrajectory')
        disp(' ');
        disp(struct(ocgTrj.octrajectory));
        disp('OdeNum:')
        disp(' ');
        disp(ocgTrj.odenum);
    else
        disp('empty ocgtrajectory');
        disp(' ');
    end
end