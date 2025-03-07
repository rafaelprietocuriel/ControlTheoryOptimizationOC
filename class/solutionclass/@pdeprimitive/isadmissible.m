function [b,adminfo]=isadmissible(pdePrim,ppdeObj,varargin)
%
% ISADMISSIBLE test if equilibrium/limit set is admissible
%
% ISADMISSIBLE(PDEPRIM,PPDEOBJ) tests if the dynprimitive object PDEPRIM or
% cell of dynprimitives is admissible in terms of the underlying optimal
% control problem PPDEOBJ. If PDEPRIM is an equilibrium it is also checked if
% it satisfies the zero condition of the dynamics.
%
% ISADMISSIBLE(PDEPRIM,PPDEOBJ,OPT) with OPT providing the tolerance for
% admissibility.
%
% ISADMISSIBLE(PDEPRIM,PPDEOBJ,OPT,USERFUNC) with USERFUNC providing a
% user defined function (main function name, which is then build from the
% model name + USERFUNC)  
%
% B = ISADMISSIBLE(PDEPRIM,PPDEOBJ,OPT) if the criteria are satisfied B=1
% otherwise B=0.
%
% [B INFO]= ISADMISSIBLE(PDEPRIM,PPDEOBJ,OPT) INFO is a structure with the
% examined values. 

opt=[];
userfunc=[];
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

if isempty(ppdeObj) || isempty(pdePrim)
    b=[];
    adminfo=[];
    return
end
AdmissibleTolerance=getocoptions(opt,'GENERAL','AdmissibleTolerance');
ZeroDeviationTolerance=getocoptions(opt,'GENERAL','ZeroDeviationTolerance');
ImaginaryTolerance=getocoptions(opt,'GENERAL','ImaginaryTolerance');
LegendreClebsch=strcmpi(getocoptions(opt,'GENERAL','LegendreClebsch'),'on');

switch modeltype(ppdeObj)
    case {'ppdemodel'}
        constraintval=admissible(ppdeObj,pdePrim);
        %zeroval=canonicalsystem(ppdeObj,pdePrim);
        zeroval=equilibriumequation(ppdeObj,pdePrim);
end
imagval=imag(pdePrim);
if ~isempty(userfunc)% && exist(userfunc,'file')
    userval=feval(ppdeObj,userfunc,independentvar(pdePrim),dependentvar(pdePrim),parametervalue(ppdeObj),arcargument(pdePrim));
else
    userval=[];
end

if LegendreClebsch
    eigval=eig(d2hamiltoniandu2(ppdeObj,pdePrim));
    leclc=all(eigval<0);
end
% b=all(constraintval(:)>=-AdmissibleTolerance) && all(abs(zeroval(:))<=ZeroDeviationTolerance*norm(pdePrim)) ...
%     && all(userval(:)>=-AdmissibleTolerance) && all(abs(imagval(:))<=ImaginaryTolerance);
b=all(constraintval(:)>=-AdmissibleTolerance) && all(abs(zeroval(:))<=ZeroDeviationTolerance) ...
    && all(userval(:)>=-AdmissibleTolerance) && all(abs(imagval(:))<=ImaginaryTolerance);

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
