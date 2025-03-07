function [varargout]=contlc(ocObj,dynPrim,ap,r0,ntst,ncol,opt,varargin)
%
% CONTLC continues a limit cycle using MATCONT
%
% CONTLC(OCOBJ,OCH,AP) usually the continuation of a limit cycle starts
% from an equilibrium at a Hopf bifurcation, detected by a previous MATCONT
% equilibrium continuation. Then OCH is the corresponding equilibrium
% stored as a dynprimitive object and AP the active parameter index.
%
% CONTLC(OCOBJ,OCH,AP,L0,NTST,NCOL) then L0 contains the value of the
% initial amplitude and NTST and NCOL are the number of mesh and
% collocation points  to be used for the discretization.
%
% CONTLC(OCOBJ,X,AP,V,S) to start a continuation from a previous limit
% cycle continuation the data X, thangent data V and information structure
% S has to be provided.  
%
% CONTLC(OCOBJ,X,AP,V,S,L0,NTST,NCOL)
%
% CONTLC(OCOBJ,X,AP,...,OPT) if the last argument is an ocmat option
% structure the options for MATCONT for the MATCONT continuation.
%
% [X V S H F] = CONTLC(...) the arguments X, V, S are the same as for the
% equilibrium continuation (CONTEP)
% F contains different information, depending on the continuation run. For
% noncycle-related continuations, the f-vector just contains the
% eigenvalues, if asked for. For limit cycle continuations, it begins with
% the mesh points of the time- discretization, followed by, if they were
% asked for, the PRC- and dPRC-values in all points of the periodic orbit
% (cf. MATCONT manual). Then, if required, follow the multipliers and
% monodromy matrix, if instead of the original "limitcycle.m" file the
% ocmat version is used (this is the default behavior).
global cds lds eds

if isequilibrium(dynPrim)
    if nargin<3
        r0=[];
        ntst=[];
        ncol=[];
        opt=[];
    end
    if nargin<4
        r0=[];
        ncol=[];
        opt=[];
    end
    if nargin<5
        r0=[];
        opt=[];
    end
    if nargin<6
        opt=[];
    end
    if isempty(r0)
        r0=1e-5;
    end
    if isempty(ntst)
        ntst=20;
    end
    if isempty(ncol)
        ncol=4;
    end
    if isempty(opt)
        opt=defaultocoptions;
    end
    if ischar(ap)
        ap=parameterindex(ocObj,ap);
    end

    eds.P0=parametervalue(ocObj);
    eds.ActiveParams=ap;
    [xlc0,vlc0]=ocmatinit_H_LC(ocObj,dynPrim,ap,r0,ntst,ncol,varargin{:});
else
    if nargin<3
        opt=[];
    else
        opt=r0;
    end
    if isempty(opt)
        opt=defaultocoptions;
    end
    [xlc0,vlc0]=ocmatinit_LC_LC(ocObj,dynPrim,ap,opt,varargin{:});
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
matcontpath=getocmatpath('cl_matcont');
if ~isempty(matcontpath)
    mydir=cd;
    cd(fullfile(matcontpath,'LimitCycle'))
    [outArgument.x,outArgument.v,outArgument.s,outArgument.h,outArgument.f]=cont(@limitcycle,xlc0,vlc0,opt.MATCONT);
    cd(mydir)
else
    outArgument.x=[];
    outArgument.v=[];
    outArgument.s=[];
    outArgument.h=[];
    outArgument.f=[];
end

globalVariable.General=cds;
globalVariable.Specific=lds;

varargout=struct2cell(outArgument);
%savematcontoutput(ocObj,globalVariable,outArgument)