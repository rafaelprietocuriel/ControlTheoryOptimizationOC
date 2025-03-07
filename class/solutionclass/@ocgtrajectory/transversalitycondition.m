function out=transversalitycondition(ocgTrj,varargin)
%
% TRANSVERSALITYCONDITION

if nargin==1
    ocObj=stdocmodel(modelname(ocgTrj),[],[],0);
end

if nargin>=2
    ocObj=varargin{1};
end
par=modelparameter(ocgTrj);
ocObj=changeparametervalue(ocObj,par);

arcarg=arcargument(ocgTrj);
indepvar=time(ocgTrj);
depvar=dependentvar(ocgTrj);


out=feval(ocObj,'TransversalityBC',indepvar(end),depvar{end}(:,end),par,arcarg);
