function [b adminfo]=isadmissible(ocCuv,ocObj,varargin)
%
% ISADMISSIBLE test if occurve object set is admissible
%
% ISADMISSIBLE(OCCURVE,OCOBJ) tests if the dynprimitive object DYNPRIM or
% cell of dynprimitives is admissible in terms of the underlying optimal
% control problem OCOBJ. If DYNPRIM is an equilibrium it is also checked if
% it satisfies the zero condition of the dynamics.
%
% ISADMISSIBLE(OCOBJ,DYNPRIM,OPT) with OPT providing the tolerance for
% admissibility.
%
% ISADMISSIBLE(OCOBJ,DYNPRIM,OPT,USERFUNC) with USERFUNC providing a
% user defined function (main function name, which is then build from the
% model name + USERFUNC)  
%
% B = ISADMISSIBLE(OCOBJ,DYNPRIMJ,OPT) if the criteria are satisfied B=1
% otherwise B=0.
%
% [B INFO]= ISADMISSIBLE(OCOBJ,DYNPRIMJ,OPT) INFO is a structure with the
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

if isempty(ocObj)
    return
end
AdmissibleTolerance=getocoptions(opt,'GENERAL','AdmissibleTolerance');
ImaginaryTolerance=getocoptions(opt,'GENERAL','ImaginaryTolerance');
LegendreClebsch=strcmpi(getocoptions(opt,'GENERAL','LegendreClebsch'),'on');

constraintval=admissible(ocObj,ocCuv);
imagval=imag(ocCuv);
if ~isempty(userfunc)
    userval=feval(ocObj,userfunc,independentvar(ocCuv),dependentvar(ocCuv),parametervalue(ocObj),arcargument(ocCuv));
else
    userval=[];
end

if LegendreClebsch
    depvar=dependentvar(ocCuv);
    arcarg=arcargument(ocCuv);
    numpts=size(depvar,2);
    leclc=zeros(1,numpts);
    for ii=1:numpts
        sol.dependentvar=depvar(:,ii);
        sol.arcarg=arcarg;
        sol.independentvar=0;
        sol.arcposition=[1 1]';
        eigval=eig(d2hamiltoniandu2(ocObj,sol));
        leclc(ii)=all(eigval<0);
    end
end
b=all(constraintval>=-AdmissibleTolerance) & (userval>=-AdmissibleTolerance) ...
    & all(abs(imagval)<=ImaginaryTolerance);

if LegendreClebsch
    b=b & leclc;
    adminfo.LegendreClebsch=eigval;
else
    adminfo.LegendreClebsch=[];
end
adminfo.ConstraintValue=constraintval;
adminfo.UserValue=userval;
adminfo.ImaginaryValue=imagval;
