function solp=partialarc(sol,partarccoord,timeshift)
%
% PARTIALARC returns a part of a trajectory
%
% PARTIALARC(SOL,PARTARCCOORD) the arc of the solution structure SOl is
% returned for the indices specified by PARTARCCOORD.
%
% SOLP=PARTIALARC(SOL,PARTARCCOORD) SOLP is a solution structure, with
% adapted fields 'x', 'y', 'arcarg', 'arcinterval', 'arcposition'.

if nargin==2
    timeshift=1;
end
solp=[];
arcposition=find(diff(sol.x)==0);
leftarcindex=[1 arcposition+1];
%rightarcindex=[arcposition numel(sol.x)];
absx=transform2absx(sol);

ff=find(leftarcindex-partarccoord(1)<=0,1,'first');
if leftarcindex(ff)-partarccoord(end)>0
    return
end

solp.x=sol.x(partarccoord);
solp.y=sol.y(:,partarccoord);
solp.arcarg=sol.arcarg(ff);
solp.arcinterval=absx([partarccoord(1) partarccoord(end)]);
solp.arcposition=[1 numel(partarccoord)];
 if timeshift
     solp.x=(solp.x-solp.x(1));
     solp.x=solp.x/solp.x(end);
     solp.arcinterval=solp.arcinterval-solp.arcinterval(1);
 end
function absx=transform2absx(sol)
% transform x (relative/normalized time) into absolute time
arcposition=find(diff(sol.x)==0);
leftarcindex=[1 arcposition+1];
rightarcindex=[arcposition numel(sol.x)];

absx=sol.x;
for ii=1:numel(leftarcindex)
    xint=sol.x(leftarcindex(ii):rightarcindex(ii));
    absx(leftarcindex(ii):rightarcindex(ii))=(xint-xint(1))*(sol.arcinterval(ii+1)-sol.arcinterval(ii))+sol.arcinterval(ii);
end