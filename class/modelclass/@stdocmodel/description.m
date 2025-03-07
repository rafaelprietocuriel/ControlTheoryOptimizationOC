function varargout=description(ocObj)

if isempty(ocObj)
    if nargout==1
        varargout{1}='';
    end
    return
end
if nargout==1
    varargout{1}=description(ocObj.Model);
    return
end

initfn=initfilename(ocObj);
fprintf('Initialization file:\n\t')
fprintf(['<a href="matlab:edit(' inputname(1) ')">' initfn '</a>'])
fprintf('\n\n')
fprintf('Name of the model:\n\t%s\n\n',modelname(ocObj))
description(ocObj.Model);
