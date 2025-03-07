function display(ocTrj)

if isequal(get(0,'FormatSpacing'),'compact')
    disp([inputname(1) ' =']);
    if ~isempty(ocTrj)
        disp('ocmatclass: hybridoctrajectory')
        disp(struct(ocTrj));
    else
        disp('empty hybridoctrajectory');
    end
else
    disp(' ')
    disp([inputname(1) ' =']);
    disp(' ')
    if ~isempty(ocTrj)
        disp('ocmatclass: hybridoctrajectory')
        disp(' ');
        disp(struct(ocTrj));
    else
        disp('empty hybridoctrajectory');
        disp(' ');
    end
end