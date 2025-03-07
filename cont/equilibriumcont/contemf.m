function varargout=contemf(ocObj,ocEP,acoord,opt,varargin)
%
% CONTEP continues an equilibrium using MATCONT
%
% CONTEP(OCOBJ,OCEP,AP) the equilibrium stored as a dynprimitve object OCEP
% is continued by changing the parameter value with index AP, starting at
% the actual parameter value of OCOBJ. The necessary MATCONT file (see
% MATCONT manual) is generated automatically as "model4matcont.m", where
% model is the name of OCOBJ.
%
% X = CONTEP(OCOBJ,OCEP,AP) the returned argument X is a matrix, where each
% column contains the values of the equilibrium along the continuation
% curve. In the last row the actual value of the continued parameter is
% written.
%
% [X V S] = CONTEP(OCOBJ,OCEP,AP) the matrix V contains the tangent vectors
% along the continuation curve and S is a structure with fields:
%   s.index  ... index of the singularity point in x, so s(1).index is
%           always equal to 1 and s(end).index is the number of computed
%           points.
%   s.label  ... label of the singularity; by convention s(1).label is ”00”
%           and s(end).label is "99".
%   s.msg    ... a string that may contain any information which is useful
%           for the user, for example the full name of the detected special
%           point. 
% For further information about the returned arguments the reader is
% referred to the MATCONT manual.
%
% [X V S] = CONTEP(OCOBJ,OCEP,AP,OPT) OPT is an occont option structure,
% where the field OPT.MATCONT provides specific options for MATCONT. To
% change MATCONT specific options use the command:
%
%   opt=setocoptions(opt,'MATCONT',OPTIONNAME,OPTIONVALUE)
%
% The most important fields for usual applications are:
%   InitStepsize ... the initial stepsize (default: 0.01)
%   MaxStepsize  ... the maximum stepsize (default: 0.1)
%   MaxNumPoints ...
%   Backward     ... boolean indicating the direction of the continuation
%               (sign of the initial tangent vector) v0 (default: 0)
%   Singularity  ... boolean indicating the presence of a singularity
%               matrix (default: 1)
% 
% for a complete list of options the reader is referred to the MATCONT
% manual
global cds eds

if nargin==3
    opt=defaultocoptions;
end
if strcmpi(getocoptions(opt,'OCCONTARG','CheckAdmissibility'),'on')
    UserInfo.name='testadmissibility';
    UserInfo.state=1;
    UserInfo.label='adm';
    opt=setocoptions(opt,'MATCONT','Userfunctions',1,'UserfunctionsInfo',UserInfo);
end

if ischar(acoord)
    acoord=parameterindex(ocObj,acoord);
end
[x0,v0]=ocmatinit_EP_EMF(ocObj,ocEP,acoord,varargin{:});
matcontpath=getocmatpath('cl_matcont');
if ~isempty(matcontpath)
    mydir=cd;
    cd(fullfile(matcontpath,'Equilibrium'))
    [outArgument.x,outArgument.v,outArgument.s,outArgument.h,outArgument.f]=cont(@equilibrium,x0,v0,opt.MATCONT);
    cd(mydir)
else
    outArgument.x=[];
    outArgument.v=[];
    outArgument.s=[];
    outArgument.h=[];
    outArgument.f=[];
end

globalVariable.General=cds;
globalVariable.Specific=eds;

varargout=struct2cell(outArgument);
%savematcontoutput(ocObj,globalVariable,outArgument)