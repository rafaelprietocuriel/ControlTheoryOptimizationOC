function disp(ocTrj)

if isequal(get(0,'FormatSpacing'),'compact')
    if ~isempty(ocTrj)
        disp('ocmatclass: hybridoctrajectory')
        disp(struct(ocTrj));
    else
        disp('empty hybridoctrajectory');
    end
else
    if ~isempty(ocTrj)
        disp('ocmatclass: hybridoctrajectory')
        disp(' ');
        disp(struct(ocTrj));
    else
        disp('empty hybridoctrajectory');
        disp(' ');
    end
end