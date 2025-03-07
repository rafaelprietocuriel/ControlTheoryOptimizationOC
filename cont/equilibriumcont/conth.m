function varargout=conth(ocObj,ocEP,ap,opt,varargin)
%
% CONTH continues a Hopf bifurcation curve using MATCONT
%

global cds hds

if nargin==3
    opt=defaultocoptions;
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

if ischar(ap)
    ap=parameterindex(ocObj,ap);
end
[x0,v0]=ocmatinit_H_H(ocObj,ocEP,ap,varargin{:});
matcontpath=getocmatpath('cl_matcont');
if ~isempty(matcontpath)
    %mydir=cd;
    %cd(fullfile(matcontpath,'Hopf'))
    [outArgument.x,outArgument.v,outArgument.s,outArgument.h,outArgument.f]=cont(@hopf,x0,v0,opt.MATCONT);
    %cd(mydir)
else
    outArgument.x=[];
    outArgument.v=[];
    outArgument.s=[];
    outArgument.h=[];
    outArgument.f=[];
end

globalVariable.General=cds;
globalVariable.Specific=hds;

varargout=struct2cell(outArgument);
%savematcontoutput(ocObj,globalVariable,outArgument)