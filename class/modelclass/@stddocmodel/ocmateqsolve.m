function [x,fval,exitflag,info,jacob]=ocmateqsolve(docObj,zerofun,t,x0,arcarg,opt,varargin)
%
%
x=[];
fval=[];
exitflag=[];
info=[];
jacob=[];
par=parametervalue(docObj);
equationsolver=getocoptions(opt,'GENERAL','EquationSolver');

switch equationsolver
    case 'fsolve'
        % test if optimization toolbox is registered
        if isempty(ver('optim'))
            ocmaterror('Optimization toolbox is not registered. Please provide another nonlinear eqaution solver.\n')
            return
        end
        [x,fval,exitflag,info,jacob]=fsolve(@zerofunction,x0,opt.EQ,t,par,arcarg,str2func([modelname(docObj) zerofun]),varargin{:});
        info.solver='fsolve';
end

function out=zerofunction(x,t,par,arcarg,zerofun,varargin)

if isempty(t)
    out=zerofun(x,par,arcarg,varargin{:});
else
    out=zerofun(t,x,par,arcarg,varargin{:});
end