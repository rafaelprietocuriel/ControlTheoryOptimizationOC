function ocMultiPath=deval(ocMultiPath,xint)
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

sol=ocmultipath2nativematlab(ocMultiPath);
ocMultiPathStruct=struct(ocMultiPath);
ocMultiPathStruct.solver='';
switch sol.solver
    case 'bvp4c'
        [Sxint,Spxint]=deval_bvp4c(sol.x,sol.y,xint,sol.yp);
        ocMultiPathStruct.y=Sxint;
        ocMultiPathStruct.x=xint;
        ocMultiPathStruct.solverinfo.yp=Spxint;
    case 'bvp5c'
        [Sxint,Spxint]=deval_bvp5c(sol.x,sol.y,xint,sol.yp,sol.ypmid);
        ocMultiPathStruct.y=Sxint;
        ocMultiPathStruct.x=xint;
        ocMultiPathStruct.solverinfo.yp=Spxint;
        ocMultiPathStruct.solverinfo=rmfield(ocMultiPathStruct.solverinfo,'ypmid');
    case 'bvp6c'
        [Sxint,Spxint]=deval_bvp6c(sol.x,sol.y,xint,sol.yp,sol.ypmid);
        ocMultiPathStruct.y=Sxint;
        ocMultiPathStruct.x=xint;
        ocMultiPathStruct.solverinfo.yp=Spxint;
        ocMultiPathStruct.solverinfo=rmfield(ocMultiPathStruct.solverinfo,'ypmid');
    otherwise
        Sxint=deval(sol.x,sol.y);
        ocMultiPathStruct.y=Sxint;
        ocMultiPathStruct.x=xint;
end
arcposition=find(diff(xint)==0);
leftarcindex=[1 arcposition+1];
rightarcindex=[arcposition numel(xint)];
if numel(arcposition)+1~=numel(arcargument(ocMultiPath))
    ocmaterror('Interpolation grid is not concistent with arc structure\n')
end
ocMultiPathStruct.arcposition=[leftarcindex;rightarcindex];
ocMultiPath=octrajectory(ocMultiPathStruct);
