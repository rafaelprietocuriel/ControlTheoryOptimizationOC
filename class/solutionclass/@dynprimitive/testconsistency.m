function b=testconsistency(dynPrim,varargin)

b=1;
opt=[];
ocObj=[];
if nargin>=2
    ocObj=varargin{1};
end
if nargin>=3
    opt=varargin{2};
end
if isempty(opt)
    opt=defaultocoptions;
end
if isempty(dynPrim)
    return
end
tol=getocoptions(opt,'GENERAL','ZeroDeviationTolerance');
arcarg=arcargument(dynPrim);
arcargnum=numel(arcarg);
sizedependentvar=size(dependentvar(dynPrim));
sizeindependentvar=size(independentvar(dynPrim));
periodval=period(dynPrim);
arcpos=arcposition(dynPrim);
sizearcpos=size(arcpos);
arcintval=arcinterval(dynPrim);
sizearcintval=size(arcintval);
if sizeindependentvar(1)~=1
    b=0;
    ocmatmsg('Number of indpendent variables larger than one.\n')
    return
end
if sizeindependentvar(2)~=sizedependentvar(2)
    b=0;
    ocmatmsg('Size of timegrid and points are different.\n')
    return
end
if periodval==0 && prod(sizeindependentvar)~=1    
    b=0;
    ocmatmsg('Size of timegrid for an equilibrium is different from one.\n')
    return
end
if sizearcpos(1)~=2
    b=0;
    ocmatmsg('Wrong size of arc position entry.\n')
    return
end
if sizearcintval(1)~=1
    b=0;
    ocmatmsg('Wrong size of arc interval entry.\n')
    return
end
if arcargnum~=sizearcpos(2)
    b=0;
    ocmatmsg('Number of arcs and arc positions are different.\n')
    return
end
if arcargnum~=sizearcintval(2)-1
    b=0;
    ocmatmsg('Number of arcs and arc intervals are different.\n')
    return
end
if ~isempty(ocObj) 
    % test if the model name and parameter values are equal
    if ~isempty(modelname(dynPrim)) && ~strcmp(modelname(ocObj),modelname(dynPrim))
        b=0;
        ocmatmsg('Models are different.\n')
        return
    end
    if ~isempty(modelparameter(dynPrim))
        try
            if norm(parametervalue(ocObj)-modelparameter(dynPrim),inf)>tol
                b=0;
                ocmatmsg('Parameter value for model and ''%s %s'' are different.\n',dynprimitiveclass(dynPrim),inputname(1))
                return
            end
        catch
            b=0;
            ocmatmsg('Parameter value for model and ''dynprimitive'' are different.\n')
            return
        end
    end
    if sizedependentvar(:,1)~=dependentvariablenum(ocObj,arcarg);
        b=0;
        ocmatmsg('Number of equations and variables are not equal.\n')
        return
    end
end