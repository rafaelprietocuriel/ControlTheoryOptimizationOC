function b=testconsistency(ppdePrim,varargin)

b=1;
opt=[];
ppdeObj=[];
if nargin>=2
    ppdeObj=varargin{1};
end
if nargin>=3
    opt=varargin{2};
end
if isempty(opt)
    opt=defaultocoptions;
end
if isempty(ppdePrim)
    return
end
tol=getocoptions(opt,'GENERAL','ZeroDeviationTolerance');

if ~isempty(ppdeObj) 
    % test if the model name and parameter values are equal
    if ~isempty(modelname(ppdePrim)) && ~strcmp(modelname(ppdeObj),modelname(ppdePrim))
        b=0;
        ocmatmsg('Model name of ''ppdeprimitive'' (%s) and oc model (%s) are different.\n',modelname(ppdePrim),modelname(ppdeObj))
        return
    end
    if ~isempty(modelparameter(ppdePrim))
        try
            if norm(parametervalue(ppdeObj)-modelparameter(ppdePrim),inf)>tol
                b=0;
                ocmatmsg('Parameter value for model and ''ppdeprimitive %s'' are different.\n',inputname(1))
                return
            end
        catch
            b=0;
            ocmatmsg('Parameter value for model and ''ppdeprimitive'' are different.\n')
            return
        end
    end
end