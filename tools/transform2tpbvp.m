function soltp=transform2tpbvp(sol,varargin)
%
% TRANSFORM2TPBVP transforms solution structure of multipoint BVP to
% twopoint BVP 
%
%
n=[];

if nargin>=2
    n=varargin{1};
end
if isempty(n)
    n=1;
end
soltp=sol;

arcposition=find(diff(sol.x)==0);
leftarcindex=[1 arcposition+1];
rightarcindex=[arcposition numel(sol.x)];
numarc=numel(rightarcindex);
soltp.y=[];

xtotal=[];
for ii=1:numarc
    xint=unique(sol.x([leftarcindex(ii) leftarcindex(ii):n:rightarcindex(ii) rightarcindex(ii)]));
    xtotal=[xtotal (xint-xint(1))/(xint(end)-xint(1))];
end
xtotal=unique(xtotal);
xtotal((abs(diff(xtotal)))<1000*eps)=[];
numtotal=numel(xtotal);
xnew=xtotal;
rowcounter=0;
idxcounter=0;
for ii=1:numarc
    rowcounter_start=rowcounter+1;
    idxcounter_start=idxcounter+1;
    idxcounter=idxcounter+numtotal;
    xint=sol.x(leftarcindex(ii):rightarcindex(ii));
    yint=sol.y(:,leftarcindex(ii):rightarcindex(ii));
    xnew(idxcounter_start:idxcounter)=(xint(end)-xint(1))*xtotal+xint(1);
    rowcounter=rowcounter+size(yint,1);
    arcposition(1:2,ii)=[rowcounter_start;rowcounter];
end
% evaluate solution on a grid that has an equal number of points at each arc

switch sol.solver
    case 'bvp4c'
        ynew=deval_bvp4c(sol.x,sol.y,xnew,sol.solverinfo.yp);
        soltp.solverinfo=rmfield(soltp.solverinfo,'yp');
    case 'bvp5c'
        ynew=deval_bvp5c(sol.x,sol.y,xnew,sol.solverinfo.yp,sol.solverinfo.ymid);
        soltp.solverinfo=rmfield(soltp.solverinfo,{'yp','ymid'});
    case 'bvp6c'
        ynew=deval_bvp6c(sol.x,sol.y,xnew,sol.solverinfo.yp,sol.solverinfo.ypmid);
        soltp.solverinfo=rmfield(soltp.solverinfo,{'yp','ypmid'});
    otherwise
        ynew=zeros(size(sol.y,1),numtotal*numarc);
        idxcounter=0;
        for ii=1:numarc
            idxcounter_start=idxcounter+1;
            idxcounter=idxcounter+numtotal;
            xint=sol.x(leftarcindex(ii):rightarcindex(ii));
            xnewint=xnew(idxcounter_start:idxcounter);
            yint=sol.y(:,leftarcindex(ii):rightarcindex(ii));
            ynew(:,idxcounter_start:idxcounter)=interp1(xint,yint.',xnewint).';
        end
end

idxcounter=0;
for ii=1:numarc
    idxcounter_start=idxcounter+1;
    idxcounter=idxcounter+numtotal;
    yint=ynew(:,idxcounter_start:idxcounter);
    soltp.y=[soltp.y;yint];
end
soltp.x=xtotal;
% arcposition now contains the corresponding rows of each arc
soltp.arcposition=arcposition;
soltp.solver='';
soltp.solverinfo.coeff=[];
soltp.solverinfo.tangent=[];
soltp.solverinfo.multiarccalc=false;
