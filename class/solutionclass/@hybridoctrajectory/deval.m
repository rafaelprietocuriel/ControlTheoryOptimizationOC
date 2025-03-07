function hocTrj=deval(hocTrj,xint)
%
% DEVAL evaluate an octrajectory
%
% DEVAL(hocTrj,XINT) evaluates the instance hocTrj of an octrajectory at the
% grid XINT. XINT has to consist at least of the initial and endtime of
% hocTrj. If hocTrj consists of multiple arcs XINIT has to contain also the
% grid points of the switching points between the different arcs.
%
% OCTRJN=DEVAL(hocTrj,XINT) OCTRJN is an instance of octrajectory evaluated
% at XINT. The 'solver' field is set to empty.

sol=hybridoctrajectory2nativematlab(hocTrj);
ocTrjStruct=struct(hocTrj);
ocTrjStruct.solver='';
switch sol.solver
    case 'bvp4c'
        [Sxint,Spxint]=deval_bvp4c(sol.x,sol.y,xint,sol.yp);
        ocTrjStruct.y=Sxint;
        ocTrjStruct.x=xint;
        ocTrjStruct.solverinfo.yp=Spxint;
    case 'bvp5c'
        [Sxint,Spxint]=deval_bvp5c(sol.x,sol.y,xint,sol.yp,sol.ypmid);
        ocTrjStruct.y=Sxint;
        ocTrjStruct.x=xint;
        ocTrjStruct.solverinfo.yp=Spxint;
        ocTrjStruct.solverinfo=rmfield(ocTrjStruct.solverinfo,'ypmid');
    case 'bvp6c'
        [Sxint,Spxint]=deval_bvp6c(sol.x,sol.y,xint,sol.yp,sol.ypmid);
        ocTrjStruct.y=Sxint;
        ocTrjStruct.x=xint;
        ocTrjStruct.solverinfo.yp=Spxint;
        ocTrjStruct.solverinfo=rmfield(ocTrjStruct.solverinfo,'ypmid');
    otherwise
        %Sxint=deval(xint,sol);
        Sxint=interp1(sol.x,sol.y.',xint).';
        ocTrjStruct.y=Sxint;
        ocTrjStruct.x=xint;
end
arcposition=find(diff(xint)==0);
leftarcindex=[1 arcposition+1]+1;
rightarcindex=[arcposition numel(xint)]+1;
ocTrjStruct.x=[sol.x0 ocTrjStruct.x sol.xT];
ocTrjStruct.y=[sol.y0 ocTrjStruct.y sol.yT];
ocTrjStruct.arcposition=[leftarcindex;rightarcindex];
hocTrj=hybridoctrajectory(ocTrjStruct);
