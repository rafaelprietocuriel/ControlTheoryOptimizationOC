function out=impulseobjective(ocObj,solObj,varargin)
%

if isempty(ocObj)
    out=[];
    return
end
if nargin==1
    solObj=[];
end
par=parametervalue(ocObj);
if ishybridoctrajectory(solObj)
    jumparg=jumpargument(solObj);
    arcarg=arcargument(solObj);
    arcintv=arcinterval(solObj);
    jdepvar=jumpdependentvar(solObj);
    depvar=dependentvar(solObj);
elseif isnumeric(solObj)
    arcarg=0;
    jumparg=0;
    indepvar=0;
    depvar=[solObj solObj];
end

if isempty(solObj)
    % return symbolic expressions
    % if repflag=1 control values/Lagrange multipliers are replaced by the
    % corresponding symbolic expressions for the optimal values
    repflag=[];
    arcarg=[];
    if nargin>=3
        jumparg=varargin{1};
    end
    if nargin>=4
        repflag=varargin{2};
    end
    if isempty(arcarg)
        jumparg=0;
    end
    if isempty(repflag)
        repflag=0;
    end
    out=feval(ocObj,'SymbolicImpulseObjectiveFunction',jumparg,repflag);
else
    out=0;
    for ii=1:length(jumparg)
        if jumparg(ii)
            out=out+feval(ocObj,'DiscountedImpulseObjectiveFunction',arcintv(ii),jdepvar(:,2*ii-1:2*ii),par,[],jumparg(ii),varargin{:});
        end
    end
end