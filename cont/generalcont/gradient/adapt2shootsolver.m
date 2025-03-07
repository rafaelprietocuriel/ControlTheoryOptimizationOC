function [coeff,tangent,extremal]=adapt2shootsolver(solinit)

global OCSHOOTCONT

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

%OCSHOOTCONT.HE.freeparameternum=numel(solinit.parameters);%additional/continuation parameters

coeff=solinit.parameters;
coeff=coeff(:);
extremal=solinit.extremal;

OCSHOOTCONT.HE.numdvariables=numel(solinit.parameters);

OCSHOOTCONT.TargetValueNum=1;

OCSHOOTCONT.HE.unitvector=[];
OCSHOOTCONT.HE.unitvector(OCSHOOTCONT.HE.numdvariables,1)=1;

switch OCSHOOTCONT.problem_func
    case 'extremalgrad2ep'
        OCSHOOTCONT.infinitetimeendconditions=1;
    otherwise
        OCSHOOTCONT.infinitetimeendconditions=0;
        
end
%OCSHOOTCONT.HE.jacobian=zeros(OCSHOOTCONT.HE.numdvariables);
OCSHOOTCONT.TIMEMESH.num=length(solinit.extremal(1).t);

