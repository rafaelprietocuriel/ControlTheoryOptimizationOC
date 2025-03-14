function varargout=contlpc(ocObj,ocLC,ap,opt,varargin)
%
% CONTLP continues a fold using MATCONT
%
% CONTLP(OCOBJ,OCEP,AP) the fold stored as a dynprimitve object OCEP
% is continued along the parameter values with index AP, starting at
% the actual parameter value of OCOBJ. The necessary MATCONT file (see
% MATCONT manual) is generated automatically as "model4matcont.m", where
% model is the name of OCOBJ.
%
% X = CONTLP(OCOBJ,OCEP,AP) the returned argument X is a matrix, where each
% column contains the values of the equilibrium along the continuation
% curve. In the last row the actual value of the continued parameter is
% written.
%
% [X V S] = CONTLP(OCOBJ,OCEP,AP) the matrix V contains the tangent vectors
% along the continuation curve and S is a structure with fields:
%   s.index  ... index of the singularity point in x, so s(1).index is
%           always equal to 1 and s(end).index is the number of computed
%           points.
%   s.label  ... label of the singularity; by convention s(1).label is �00�
%           and s(end).label is "99".
%   s.msg    ... a string that may contain any information which is useful
%           for the user, for example the full name of the detected special
%           point. 
% For further information about the returned arguments the reader is
% referred to the MATCONT manual.
%
% [X V S] = CONTLP(OCOBJ,OCEP,AP,OPT) OPT is an occont option structure,
% where the field OPT.MATCONT provides specific options for MATCONT. To
% change MATCONT specific options use the command:
%
%   opt=setocoptions(opt,'MATCONT',OPTIONNAME,OPTIONVALUE)
%
% The most important fields fpr usual applications are:
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
global cds lcpds

if nargin==3
    opt=defaultocoptions;
end
if ischar(ap)
    ap=parameterindex(ocObj,ap);
end
if getocoptions(opt,'OCCONTARG','CheckAdmissibility')
    if isempty(opt.MATCONT.Userfunctions)
        UserInfo.name='testadmissibility';
        UserInfo.state=1;
        UserInfo.label='adm';
        opt=setocoptions(opt,'MATCONT','Userfunctions',1,'UserfunctionsInfo',UserInfo);
    else
        l=length(opt.MATCONT.Userfunctions);
        UserInfo(l+1)=opt.MATCONT.UserfunctionsInfo;
        UserInfo(1).name='testadmissibility';
        UserInfo(1).state=1;
        UserInfo(1).label='adm';
        userfunc=[1 opt.MATCONT.Userfunctions];
        userfuncstop=opt.MATCONT.UserfunctionsStop;
        opt=setocoptions(opt,'MATCONT','UserfunctionsStop',userfuncstop,'Userfunctions',userfunc,'UserfunctionsInfo',UserInfo);
    end
end
[x0,v0]=ocmatinit_LPC_LPC(ocObj,ocLC,ap,varargin{:});
matcontpath=getocmatpath('cl_matcont');
if ~isempty(matcontpath)
%     mydir=cd;
%     cd(fullfile(matcontpath,'LimitPointCycle'))
    [outArgument.x,outArgument.v,outArgument.s,outArgument.h,outArgument.f]=cont(@limitpointcycle,x0,v0,opt.MATCONT);
%    cd(mydir)
else
    outArgument.x=[];
    outArgument.v=[];
    outArgument.s=[];
    outArgument.h=[];
    outArgument.f=[];
end

globalVariable.General=cds;
globalVariable.Specific=lcpds;

varargout=struct2cell(outArgument);
