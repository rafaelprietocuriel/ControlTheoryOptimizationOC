function [x,fval,exitflag,info,jacob]=ocmateqsolve(ocObj,zerofun,t,x0,arcarg,opt,varargin)
%
%
x=[];
fval=[];
exitflag=[];
info=[];
jacob=[];
par=parametervalue(ocObj);
equationsolver=getocoptions(opt,'GENERAL','EquationSolver');

switch equationsolver
    case 'fsolve'
        % test if optimization toolbox is registered
        if isempty(ver('optim'))
            ocmaterror('Optimization toolbox is not registered. Please provide another nonlinear equation solver.\n')
            return
        end
        if strcmp(opt.EQ.Jacobian,'on')
            try
                [x,fval,exitflag,info,jacob]=fsolve(@zerofunctionwithjac,x0,opt.EQ,t,par,arcarg,str2func([modelname(ocObj) zerofun]),str2func([modelname(ocObj) varargin{1}]),varargin{2:end});
                info.solver='fsolve';
            end
        else
            varargin(1)=[];
            try
                [x,fval,exitflag,info,jacob]=fsolve(@zerofunction,x0,opt.EQ,t,par,arcarg,str2func([modelname(ocObj) zerofun]),varargin{:});
                info.solver='fsolve';
            end
        end

    case 'mynewtsolve'

        if strcmp(opt.EQ.Jacobian,'on')
            jacfun=str2func([modelname(ocObj) varargin{1}]);
        else
            jacfun=[];
        end
        [x,fval,exitflag,info]=mynewtsolve(x0,str2func([modelname(ocObj) zerofun]),jacfun,opt,par,arcarg);
        info.solver='fsolve';
end

function out=zerofunction(x,t,par,arcarg,zerofun,varargin)

try
    out=zerofun(x,par,arcarg,varargin{:});
catch
    out=zerofun(t,x,par,arcarg,varargin{:});
end

function [fx,J]=zerofunctionwithjac(x,t,par,arcarg,zerofun,jacfun,varargin)

try
    fx=zerofun(x,par,arcarg,varargin{:});
    J=jacfun(x,par,arcarg,varargin{:});
catch
    fx=zerofun(t,x,par,arcarg,varargin{:});
    J=jacfun(t,x,par,arcarg,varargin{:});
end