function display(ocTrj)

if isequal(get(0,'FormatSpacing'),'compact')
    disp([inputname(1) ' =']);
    if ~isempty(ocTrj)
        disp('ocmatclass: octrajectory')
        disp(struct(ocTrj));
    else
        disp('empty octrajectory');
    end
else
    disp(' ')
    disp([inputname(1) ' =']);
    disp(' ')
    if ~isempty(ocTrj)
        disp('ocmatclass: octrajectory')
        disp(' ');
        disp(struct(ocTrj));
    else
        disp('empty octrajectory');
        disp(' ');
    end
end