function ocTrj=cat(dim,varargin)

if nargin==1
    ocTrj=octrajectory();
    return
end

switch dim
    case 1
        if isocasymptotic(varargin{1})
            ocTrj=octrajectory(varargin{1});
            ocTrj.timehorizon=inf;
        else
            ocTrj=varargin{1};
        end
        for ii=2:nargin-1
            if ocTrj.arcarg(end)~=varargin{ii}.arcarg(1)
                ocTrj=addarc(ocTrj,varargin{ii});
            else
                ocTrj=extendarc(ocTrj,varargin{ii});
            end
        end
        ocTrj.userinfo=[];
        ocTrj.violationinfo=[];
        ocTrj.solverinfo=[];
        ocTrj.solver=[];
    case 2
end

function ocTrj=addarc(ocTrj,ocTrj2)

ocTrj.arcarg=[ocTrj.arcarg ocTrj2.arcarg];
ocTrj.arcposition=[ocTrj.arcposition ocTrj.arcposition(2,end)+ocTrj2.arcposition];
if ocTrj.arcinterval(end)~=ocTrj2.arcinterval(1)
    ocmaterror('Concatenation failed. octrajectories are incompatible.')
end
ocTrj.arcinterval=[ocTrj.arcinterval ocTrj2.arcinterval(2:end)];
ocTrj.y=[ocTrj.y ocTrj2.y];
ocTrj.x=[ocTrj.x ocTrj.x(end)+ocTrj2.x];

function ocTrj=extendarc(ocTrj,ocTrj2)

if arcnum(ocTrj2)>1
    ocmaterror('Concatenation failed. octrajectories are incompatible.')
end
t=absolutetime(ocTrj);
t2=absolutetime(ocTrj2);
ocTrj.arcposition(2,end)=ocTrj.arcposition(2,end)+ocTrj2.arcposition(2,end)-1;
if ocTrj.arcinterval(end)~=ocTrj2.arcinterval(1)
    ocmaterror('Concatenation failed. octrajectories are incompatible.')
end
ocTrj.arcinterval(end)=ocTrj2.arcinterval(end);
ocTrj.y=[ocTrj.y ocTrj2.y(:,2:end)];
t=[t t2(2:end)];
x=[];
for ii=1:arcnum(ocTrj)
    tint=t(ocTrj.arcposition(1,ii):ocTrj.arcposition(2,ii));
    x=[x (tint-tint(1))/(tint(end)-tint(1))+ii-1];
end
ocTrj.x=x;