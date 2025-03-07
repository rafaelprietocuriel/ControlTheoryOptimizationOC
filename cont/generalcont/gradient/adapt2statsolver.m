function [coeff,tangent]=adapt2statsolver(solinit)

global OCSTATCONT
% Validate arguments
if ~isfield(solinit,'coeff')
    msg=sprintf('The field ''coeff'' not present in SOLINIT.');
    ocmaterror('MATLAB:bvpcont:NoXInSolinit','%s',msg);
end

if isfield(solinit,'solverinfo') && isfield(solinit.solverinfo,'tangent')
    tangent=solinit.solverinfo.tangent;
else
    tangent=[];
end

%OCSTATCONT.HE.freeparameternum=numel(solinit.parameters);%additional/continuation parameters

coeff=solinit.coeff;

OCSTATCONT.HE.numdvariables=numel(coeff);

OCSTATCONT.TargetValueNum=1;

OCSTATCONT.HE.unitvector=[];
OCSTATCONT.HE.unitvector(OCSTATCONT.HE.numdvariables,1)=1;
if ~OCSTATCONT.OPTIONS.continuationmethod
    coeff(end)=[];
end
% 
% warnoffId={'MATLAB:singularMatrix','MATLAB:nearlySingularMatrix'};
% for ii=1:length(warnoffId)
%     warnstat(ii)=warning('query',warnoffId{ii});
%     warnoff(ii)=warnstat(ii);
%     warnoff(ii).state='off';
% end
% OCSTATCONT.warnstat=warnstat;
% OCSTATCONT.warnoff=warnoff;
% OCSTATCONT.warnoffId=warnoffId;


