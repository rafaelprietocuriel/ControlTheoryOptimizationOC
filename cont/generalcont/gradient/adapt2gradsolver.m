function [coeff,tangent,extremal]=adapt2gradsolver(solinit)

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

OCGRADCONT.linestepwidth=OCGRADCONT.OPTIONS.initlinestepwidth;

OCGRADCONT.HE.unitvector=[];
OCGRADCONT.HE.unitvector(OCGRADCONT.HE.numdvariables,1)=1;

OCGRADCONT.Jacobian4Equation=zeros(length(coeff)-1,length(coeff));

switch OCGRADCONT.problem_func
    case 'extremalgrad2ep'
        OCGRADCONT.infinitetimeendconditions=1;
    otherwise
        OCGRADCONT.infinitetimeendconditions=0;
        
end
%OCGRADCONT.HE.jacobian=zeros(OCGRADCONT.HE.numdvariables);
OCGRADCONT.TIMEMESH.num=length(solinit.extremal(1).t);
if ~OCGRADCONT.cconstraint
    OCGRADCONT.OPTIONS.directioncorrectiontype='';
    OCGRADCONT.OPTIONS.projection=0;
end
if any(strcmp(OCGRADCONT.OPTIONS.gradientmappingmethod,{'explicit','nesterov'}))
    %OCGRADCONT.OPTIONS.initlinestepwidth=min(OCGRADCONT.OPTIONS.initlinestepwidth,1/OCGRADCONT.OPTIONS.gradientmappinggamma);
    OCGRADCONT.Nesterov.QuadraticMatrix=OCGRADCONT.OPTIONS.gradientmappinggamma*speye(OCGRADCONT.control_num.concentrated*OCGRADCONT.TIMEMESH.num);
end



