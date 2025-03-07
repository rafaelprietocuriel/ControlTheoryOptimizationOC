function varargout=contbranch(ocObj,ocEP,s,h,opt,varargin)
%
% CONTH continues a Hopf bifurcation curve using MATCONT
%

global cds eds


if nargin==4
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
[x0,v0]=ocmatinit_BP_EP(ocObj,ocEP,s,h,varargin{:});
matcontpath=getocmatpath('cl_matcont');
if ~isempty(matcontpath)
    mydir=cd;
    cd(fullfile(matcontpath,'BranchPoint'))
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