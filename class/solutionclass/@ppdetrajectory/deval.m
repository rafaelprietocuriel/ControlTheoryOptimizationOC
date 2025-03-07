function ppdeTrj=deval(ppdeTrj,xint)
%
% DEVAL evaluate an ppdetrajectory
%
% DEVAL(OCTRJ,XINT) evaluates the instance OCTRJ of an ppdetrajectory at the
% grid XINT. XINT has to consist at least of the initial and endtime of
% OCTRJ. If OCTRJ consists of multiple arcs XINIT has to contain also the
% grid points of the switching points between the different arcs.
%
% OCTRJN=DEVAL(OCTRJ,XINT) OCTRJN is an instance of ppdetrajectory evaluated
% at XINT. The 'solver' field is set to empty.

sol=ppdetrajectory2nativematlab(ppdeTrj);
ppdeTrjStruct=struct(ppdeTrj);
switch sol.solver
    case 'bvp4c'
        [Sxint,Spxint]=deval_bvp4c(sol.x,sol.y,xint,sol.yp);
        ppdeTrjStruct.y=Sxint;
        ppdeTrjStruct.t=xint;
        ppdeTrjStruct.discretizationinfo.yp=Spxint;
    case 'bvp5c'
        [Sxint,Spxint]=deval_bvp5c(sol.x,sol.y,xint,sol.yp,sol.ypmid);
        ppdeTrjStruct.y=Sxint;
        ppdeTrjStruct.t=xint;
        ppdeTrjStruct.discretizationinfo.yp=Spxint;
        ppdeTrjStruct.discretizationinfo=rmfield(ppdeTrjStruct.discretizationinfo,'ypmid');
    case 'bvp6c'
        [Sxint,Spxint]=deval_bvp6c(sol.x,sol.y,xint,sol.yp,sol.ypmid);
        ppdeTrjStruct.y=Sxint;
        ppdeTrjStruct.t=xint;
        ppdeTrjStruct.discretizationinfo.yp=Spxint;
        ppdeTrjStruct.discretizationinfo=rmfield(ppdeTrjStruct.discretizationinfo,'ypmid');
    otherwise
        Sxint=deval(sol.x,sol.y);
        ppdeTrjStruct.y=Sxint;
        ppdeTrjStruct.t=xint;
end
ppdeTrj=ppdetrajectory(ppdeTrjStruct);
