function a=admissible(dgObj,depvar,arcarg,varargin)
%
% ADMISSIBLE returns the value of the constraint adn Lagrangian multipliers
%
a=[];
useradmissiblefunc=[];
indepvar=0;
if isempty(dgObj) || isempty(depvar)
    return
end
if nargin<=2
    arcarg=0;
end
if isstruct(depvar)
    try
        arcarg=depvar.arcarg;
        indepvar=depvar.independentvar;
        depvar=depvar.dependentvar;
    catch
        ocmaterror('If second input argument is a structure the fields ''dependentvar'' and ''arcarg'' have to exist!')
    end
elseif isdynprimitive(depvar) || isoctrajectory(depvar) || isocasymptotic(depvar) || isoccurve(depvar)
    arcarg=arcargument(depvar);
    indepvar=independentvar(depvar);
    depvar=dependentvar(depvar);
end
% return 
if nargin>=4
    indepvar=varargin{1};
end
par=parametervalue(dgObj);
a=feval(dgObj,'Admissible',indepvar,depvar,par,arcarg);
if nargin>=4
    useradmissiblefunc=varargin{1};
end
if ischar(useradmissiblefunc)
    a=[a; ...
        feval(dgObj,useradmissiblefunc,indepvar,depvar,par,arcarg)];
end