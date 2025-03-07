function ocTrj=deval(ocTrj,xint)
%
% DEVAL evaluate an octrajectory
%
% DEVAL(OCTRJ,XINT) evaluates the instance OCTRJ of an octrajectory at the
% grid XINT. XINT has to consist at least of the initial and endtime of
% OCTRJ. If OCTRJ consists of multiple arcs XINIT has to contain also the
% grid points of the switching points between the different arcs.
%
% OCTRJN=DEVAL(OCTRJ,XINT) OCTRJN is an instance of octrajectory evaluated
% at XINT. The 'solver' field is set to empty.

sol=octrajectory2nativematlab(ocTrj);
ocTrjStruct=struct(ocTrj);
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
    case 'ode45'
        Sxint=deval(sol,xint);
        ocTrjStruct.y=Sxint;
        ocTrjStruct.x=xint;
        ocTrjStruct.x=ocTrjStruct.x/ocTrjStruct.x(end);
    otherwise
        %Sxint=deval(xint,sol);
        Sxint=interp1(sol.x,sol.y.',xint).';
        ocTrjStruct.y=Sxint;
        ocTrjStruct.x=xint;
end
arcposition=find(diff(xint)==0);
leftarcindex=[1 arcposition+1];
rightarcindex=[arcposition numel(xint)];
if numel(arcposition)+1~=numel(arcargument(ocTrj))
    ocmaterror('Interpolation grid is not concistent with arc structure\n')
end
ocTrjStruct.arcposition=[leftarcindex;rightarcindex];
ocTrj=octrajectory(ocTrjStruct);
