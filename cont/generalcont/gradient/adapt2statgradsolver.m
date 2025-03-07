function [coeff,tangent,extremal]=adapt2statgradsolver(solinit)

global OCGRADCONT
% Validate arguments
if ~isfield(solinit,'parameters')
    msg=sprintf('The field ''parameters'' not present in SOLINIT.');
    ocmaterror('MATLAB:bvpcont:NoXInSolinit','%s',msg);
elseif ~isfield(solinit,'extremal')
    msg=sprintf('The field ''extremal'' not present in SOLINIT.');
    ocmaterror('MATLAB:bvpcont:NoXInSolinit','%s',msg);
end

if isfield(solinit,'solverinfo') && isfield(solinit.solverinfo,'tangent')
    tangent=solinit.solverinfo.tangent;
else
    tangent=[];
end

%OCGRADCONT.HE.freeparameternum=numel(solinit.parameters);%additional/continuation parameters

coeff=solinit.parameters;
coeff=coeff(:);
extremal=solinit.extremal;

OCGRADCONT.HE.numdvariables=numel(solinit.parameters);

OCGRADCONT.TargetValueNum=1;

OCGRADCONT.HE.unitvector=[];
OCGRADCONT.HE.unitvector(OCGRADCONT.HE.numdvariables,1)=1;
