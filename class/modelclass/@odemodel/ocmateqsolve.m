function [x,fval,exitflag,info,jacob]=ocmateqsolve(ocObj,zerofun,x0,arcarg,opt,varargin)
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
            ocmaterror('Optimization toolbox is not registered. Please provide another nonlinear eqaution solver.\n')
            return
        end
        [x,fval,exitflag,info,jacob]=fsolve(@zerofunction,x0,opt.EQ,par,arcarg,str2func([modelname(ocObj) zerofun]),varargin{:});
        info.solver='fsolve';
end

function out=zerofunction(x,par,arcarg,zerofun,varargin)

out=zerofun(x,par,arcarg,varargin{:});