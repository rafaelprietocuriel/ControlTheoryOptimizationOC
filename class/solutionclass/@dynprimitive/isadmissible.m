function [b,adminfo]=isadmissible(dynPrim,ocObj,varargin)
%
% ISADMISSIBLE test if equilibrium/limit set is admissible
%
% ISADMISSIBLE(DYNPRIM,OCOBJ) tests if the dynprimitive object DYNPRIM or
% cell of dynprimitives is admissible in terms of the underlying optimal
% control problem OCOBJ. If DYNPRIM is an equilibrium it is also checked if
% it satisfies the zero condition of the dynamics.
%
% ISADMISSIBLE(DYNPRIM,OCOBJ,OPT) with OPT providing the tolerance for
% admissibility.
%
% ISADMISSIBLE(DYNPRIM,OCOBJ,OPT,USERFUNC) with USERFUNC providing a
% user defined function (main function name, which is then build from the
% model name + USERFUNC)  
%
% B = ISADMISSIBLE(DYNPRIM,OCOBJ,OPT) if the criteria are satisfied B=1
% otherwise B=0.
%
% [B INFO]= ISADMISSIBLE(DYNPRIM,OCOBJ,OPT) INFO is a structure with the
% examined values. 

opt=[];
userfunc=[];
dhdu=[];

if nargin<2
    ocmaterror('Number of input arguments is less than two.')
end

if nargin>=3
    opt=varargin{1};
end

if nargin>=4
    userfunc=varargin{2};
end

if isempty(opt)
    opt=defaultocoptions;
end

if isempty(ocObj) || isempty(dynPrim)
    b=[];
    adminfo=[];
    return
end

AdmissibleTolerance=getocoptions(opt,'GENERAL','AdmissibleTolerance');
ZeroDeviationTolerance=getocoptions(opt,'GENERAL','ZeroDeviationTolerance');
ImaginaryTolerance=getocoptions(opt,'GENERAL','ImaginaryTolerance');
LegendreClebsch=strcmpi(getocoptions(opt,'GENERAL','LegendreClebsch'),'on');

switch modeltype(ocObj)
    case {'standardmodel','differentialgame'}
        constraintval=admissible(ocObj,dynPrim);
        %zeroval=canonicalsystem(ocObj,dynPrim);
        %zeroval=equilibriumequation(ocObj,dynPrim);
        zeroval=equilibriumequation(ocObj,dynPrim,varargin{3:end});
        %try
        %    dhdu=dhamiltoniandu(dynPrim);
        %end
    case 'odemodel'
        constraintval=[];
        %zeroval=dynamics(ocObj,dynPrim);
        zeroval=equilibriumequation(ocObj,dynPrim,varargin{3:end});
end
imagval=imag(dynPrim);
if ~isempty(userfunc)% && exist(userfunc,'file')
    userval=feval(ocObj,userfunc,independentvar(dynPrim),dependentvar(dynPrim),parametervalue(ocObj),arcargument(dynPrim));
else
    userval=[];
end

if LegendreClebsch
    eigval=eig(d2hamiltoniandu2(ocObj,dynPrim));
    leclc=all(eigval<0);
end
% b=all(constraintval(:)>=-AdmissibleTolerance) && all(abs(zeroval(:))<=ZeroDeviationTolerance*norm(dynPrim)) ...
%     && all(userval(:)>=-AdmissibleTolerance) && all(abs(imagval(:))<=ImaginaryTolerance);
b=all(constraintval(:)>=-AdmissibleTolerance) && all(abs(zeroval(:))<=ZeroDeviationTolerance) ...
    && all(userval(:)>=-AdmissibleTolerance) && all(abs(imagval(:))<=ImaginaryTolerance);

%if ~isempty(dhdu)
%    b=b&&norm(dhdu)<1e-10;
%end
if LegendreClebsch
    b=b && leclc;
    adminfo.LegendreClebsch=eigval;
else
    adminfo.LegendreClebsch=[];
end
adminfo.ConstraintValue=constraintval;
adminfo.ZeroValue=zeroval;
adminfo.UserValue=userval;
adminfo.ImaginaryValue=imagval;
