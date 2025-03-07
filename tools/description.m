function varargout=description(ocStruct)


if isempty(ocStruct)
    if nargout==1
        varargout{1}='';
    end
    return
end
if nargout==1
    varargout{1}=ocStruct.description;
    return
end
fprintf('Model description:\n')
for ii=1:length(ocStruct.description)
    fprintf('\t%s\n',ocStruct.description{ii})
end
fprintf('\n')