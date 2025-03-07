function disp(ocgTrj)

if isequal(get(0,'FormatSpacing'),'compact')
    if ~isempty(ocgTrj)
        disp('ocmatclass: ocgtrajectory')
        disp(struct(ocgTrj.octrajectory));
        disp('OdeNum:')
        disp(ocgTrj.odenum);
    else
        disp('empty ocgtrajectory');
    end
else
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