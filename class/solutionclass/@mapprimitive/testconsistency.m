function b=testconsistency(mapPrim,varargin)

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
if isempty(mapPrim)
    return
end
tol=getocoptions(opt,'GENERAL','ZeroDeviationTolerance');
arcarg=arcargument(mapPrim);
arcargnum=numel(arcarg);
sizedependentvar=size(dependentvar(mapPrim));
sizeindependentvar=size(independentvar(mapPrim));
periodval=period(mapPrim);
arcpos=arcposition(mapPrim);
sizearcpos=size(arcpos);
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
if periodval==0 && prod(sizeindependentvar)~=2   
    b=0;
    ocmatmsg('Wrong size of timesteps for a fix point is.\n')
    return
end
if sizearcpos(1)~=2
    b=0;
    ocmatmsg('Wrong size of arc position entry.\n')
    return
end
if arcargnum~=sizearcpos(2)
    b=0;
    ocmatmsg('Number of arcs and arc positions are different.\n')
    return
end
if ~isempty(ocObj) 
    % test if the model name and parameter values are equal
    if ~isempty(modelname(mapPrim)) && ~strcmp(modelname(ocObj),modelname(mapPrim))
        b=0;
        ocmatmsg('Models are different.\n')
        return
    end
    if ~isempty(modelparameter(mapPrim))
        try
            if norm(parametervalue(ocObj)-modelparameter(mapPrim),inf)>tol
                b=0;
                ocmatmsg('Parameter value for model and ''%s %s'' are different.\n',dynprimitiveclass(mapPrim),inputname(1))
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