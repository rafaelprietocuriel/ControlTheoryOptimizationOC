function b=isconsistent(ocObj,solObj)
%
% ISCONSISTENT test if model and calculation object are consitatent with
% respect to the model name and parameter values

b=1;

try
    if ~isempty(modelname(solObj)) && ~strcmp(modelname(ocObj),modelname(solObj))
        b=0;
        return
    end
    parval1=modelparameter(solObj);
    parval2=parametervalue(ocObj);
    if ~isempty(parval1)
        if numel(parval1)~=numel(parval2) || any(abs(parval1-parval2)>1e-8)
            b=0;
            return
        end
    end
end