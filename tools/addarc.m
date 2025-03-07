function sol=addarc(sol,arc,iposition)
%
% ADDARC fits an arc into a solution
%
% ADDARC(SOL,ARC) the input arguments SOL and ARC are OCMat solution
% structures. ARC is placed as a new arc in front of SOL.
if nargin==2
    iposition=1;
end

arcposition=find(diff(sol.x)==0);
leftarcindex=[1 arcposition+1];
rightarcindex=[arcposition numel(sol.x)];
numarc=numel(rightarcindex);
% time index for interior solution
if numel(iposition)==1
    iposition=[iposition iposition];
end
ff=find(leftarcindex-iposition(1)<=0,1,'last');
if ff==1
    Ltimeindex=[];
else
    Ltimeindex=[leftarcindex(ff-1):rightarcindex(ff-1)];
end
if ff<numarc
    Rtimeindex=[leftarcindex(ff+1):rightarcindex(ff+1)];
else
    Rtimeindex=[];
end
Larcindex=1:ff-1;
Rarcindex=ff+1:numarc;
ocAsymL.x=sol.x(Ltimeindex);
ocAsymL.y=sol.y(:,Ltimeindex);
ocAsymL.arcarg=sol.arcarg(Larcindex);
ocAsymL.arcposition=sol.arcposition(:,Larcindex);
ocAsymL.arcinterval=sol.arcinterval(:,Larcindex);
ocAsymR.x=sol.x(Rtimeindex);
ocAsymR.y=sol.y(:,Rtimeindex);
ocAsymR.arcarg=sol.arcarg(Rarcindex);
ocAsymR.arcposition=sol.arcposition(:,Rarcindex);
ocAsymR.arcinterval=sol.arcinterval(:,Rarcindex+1);

% absolute position with respect to sol.x
%absiposition=iposition;
if ~isempty(Ltimeindex)
    % relative position with respect to xint
    iposition=iposition-Ltimeindex(end);
end
xint=sol.x(leftarcindex(ff):rightarcindex(ff));
xintleft=xint(1:iposition(1));
xintleft=transform2unit(xintleft);
xintmiddle=arc.x;
xintmiddle=transform2unit(xintmiddle);
xintright=xint(iposition(end):end);
xintright=transform2unit(xintright);
if isempty(xintleft)
    offsettimeshift=-1;
else
    offsettimeshift=0;
end
xmiddle=[xintleft+ff-1 xintmiddle+ff xintright+ff+1]+offsettimeshift;

yint=sol.y(:,leftarcindex(ff):rightarcindex(ff));
if ~isempty(xintleft)
    yintleft=yint(:,1:iposition(1));
    %arcpositionintleft=[sol.arcposition(1,ff);absiposition(1)];
    arcargintleft=sol.arcarg(ff);
    arcintervalintleft=sol.arcinterval(ff);
else
    yintleft=[];
    %arcpositionintleft=[];
    arcargintleft=[];
    arcintervalintleft=[];
end
if ~isempty(xintright)
    yintright=yint(:,iposition(end):end);
    %arcpositionintright=[absiposition(end)+2;sol.arcposition(2,ff)+2];
    arcargintright=sol.arcarg(ff);
    arcintervalintright=sol.arcinterval(ff+1);
else
    yintright=[];
    %arcpositionintright=[];
    arcargintright=[];
    arcintervalintright=[];
end
%arcpositionmiddle=[arcpositionintleft [absiposition(1)+1;absiposition(end)+1] arcpositionintright]+offsettimeshift;
ymiddle=[yintleft arc.y yintright];
xnew=[ocAsymL.x xmiddle ocAsymR.x+2+offsettimeshift];
ynew=[ocAsymL.y ymiddle ocAsymR.y];

arcargmiddle=[arcargintleft arc.arcarg arcargintright];
arcintervalmiddle=[arcintervalintleft arc.arcinterval arcintervalintright];
sol.x=xnew;
sol.y=ynew;
sol.arcarg=[ocAsymL.arcarg arcargmiddle ocAsymR.arcarg];
arcposition=find(diff(sol.x)==0);
leftarcindex=[1 arcposition+1];
rightarcindex=[arcposition numel(sol.x)];
sol.arcposition=[leftarcindex;rightarcindex];
%sol.arcposition=[ocAsymL.arcposition arcpositionmiddle ocAsymR.arcposition];
sol.arcinterval=[ocAsymL.arcinterval arcintervalmiddle ocAsymR.arcinterval];
sol.solver='';
sol.solverinfo.parameters=[];

%--------------------------------------------------------------------------
function xnew=transform2unit(x)
xnew=[];
n=numel(x);
if n<2
    return
end
l=x(n)-x(1);
if l==0
    xnew=linspace(0,1,n);
    return
end

xnew=(x-x(1))/l;