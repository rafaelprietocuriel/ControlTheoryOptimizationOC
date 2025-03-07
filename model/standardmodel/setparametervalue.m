function ocStruct=setparametervalue(ocStruct,parval)
%
%
if isempty(ocStruct)
    return
end
parametername=fieldnames(ocStruct.parameter.variable);
if ocStruct.parameter.num~=numel(parval)
    ocmaterror('Wrong number of parameter values.')
end

parameteralgebraictermindex=retrievemodelinformation(ocStruct,'parameteralgebraictermindex');
if isempty(parameteralgebraictermindex.value)
    for ii=1:ocStruct.parameter.num
        ocStruct.parameter.variable.(parametername{ii})=parval(ii);
    end
else
    parametervaluetermindex=retrievemodelinformation(ocStruct,'parametervaluetermindex');
    for iicounter=parametervaluetermindex.value
        ocStruct.parameter.variable.(parametername{iicounter})=parval(iicounter);
    end
end