function disp(ocTrj)

if isequal(get(0,'FormatSpacing'),'compact')
    if ~isempty(ocTrj)
        disp('ocmatclass: octrajectory')
        disp(struct(ocTrj));
    else
        disp('empty octrajectory');
    end
else
    if ~isempty(ocTrj)
        disp('ocmatclass: octrajectory')
        disp(' ');
        disp(struct(ocTrj));
    else
        disp('empty octrajectory');
        disp(' ');
    end
end