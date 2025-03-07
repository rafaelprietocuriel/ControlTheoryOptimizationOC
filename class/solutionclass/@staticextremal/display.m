function display(statEx)

for ii=1:multiplicity(statEx)
    if multiplicity(statEx)>1
        disp(['Solution ' num2str(ii)])
    end
    if isequal(get(0,'FormatSpacing'),'compact')
        disp([inputname(1) ' =']);
        if ~isempty(statEx(ii))
            disp('ocmatclass: staticextremal')
            disp(struct(statEx(ii)));
        else
            disp('empty staticextremal');
        end
    else
        disp(' ')
        disp([inputname(1) ' =']);
        disp(' ')
        if ~isempty(statEx(ii))
            disp('ocmatclass: staticextremal')
            disp(' ');
            disp(struct(statEx(ii)));
        else
            disp('empty staticextremal');
            disp(' ');
        end
    end
end