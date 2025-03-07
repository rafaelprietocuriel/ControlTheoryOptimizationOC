function display(ocgMTrj)

if isequal(get(0,'FormatSpacing'),'compact')
    disp([inputname(1) ' =']);
    if ~isempty(ocgMTrj)
        disp('ocmatclass: ocgradmultipath')
        disp(['degree   : ' int2str(degree(ocgMTrj))])
        disp(struct(ocgMTrj));
    else
        disp('empty ocgradmultipath');
    end
else
    disp(' ')
    disp([inputname(1) ' =']);
    disp(' ')
    if ~isempty(ocgMTrj)
        disp('ocmatclass: ocgradmultipath')
        disp(' ');
        disp(['degree   : ' int2str(degree(ocgMTrj))])
        disp(' ');
        disp(struct(ocgMTrj));
    else
        disp('empty ocgradmultipath');
        disp(' ');
    end
end