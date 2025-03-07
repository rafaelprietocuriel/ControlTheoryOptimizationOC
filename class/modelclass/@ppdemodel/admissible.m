function a=admissible(ppdeObj,solObj,arcarg,varargin)
%
% ADMISSIBLE returns the value of the constraint adn Lagrangian multipliers
%
a=[];
%useradmissiblefunc=[];
if isempty(ppdeObj) || isempty(solObj)
    return
end
if nargin<=2
    arcarg=0;
end
if isstruct(solObj)
    try
        arcarg=solObj.arcarg;
        indepvar=solObj.independentvar;
        depvar=solObj.dependentvar;
    catch
        ocmaterror('If second input argument is a structure the fields ''dependentvar'' and ''arcarg'' have to exist!')
    end
elseif ispdeprimitive(solObj) || ispdetrajectory(solObj) || ispdeasymptotic(solObj)
    arcarg=arcargument(solObj);
    indepvar=time(ppdeObj,solObj,1);
    femdat=femdata(solObj);
    depvar=dependentvar(solObj);
end
% return 
par=parametervalue(ppdeObj);
a=feval(ppdeObj,'Admissible',indepvar,depvar,par,arcarg,femdat);
% if nargin>=4
%     useradmissiblefunc=varargin{1};
% end
% if ischar(useradmissiblefunc)
%     a=[a; ...
%         feval(ppdeObj,useradmissiblefunc,indepvar,depvar,par,arcarg)];
% end