function sol=odestruct(ocgTrj,parflag)
%
%

sol=[];
if isempty(ocgTrj)
    return
end
if nargin==1
    parflag=false;
end


sol.x=ocgTrj.octrajectory.x;
arcposition=ocgTrj.octrajectory.arcposition;

odenum=ocgTrj.odenum;
sol.y=zeros(max(odenum),length(sol.x));
for ii=1:length(odenum)
    sol.y(1:odenum(ii),arcposition(1,ii):arcposition(2,ii))=ocgTrj.octrajectory.y{ii}(1:odenum(ii),:);
end
sol.arcinterval=ocgTrj.octrajectory.arcinterval; %
sol.arcarg=ocgTrj.octrajectory.arcarg;
sol.arcposition=arcposition;
sol.parameters=[];
if parflag
    if isfield(ocgTrj.octrajectory.solverinfo,'parameters')
        sol.parameters=ocgTrj.octrajectory.solverinfo.parameters;
    end
end
sol.x0=0;
sol.solverinfo.tangent=[];
sol.solverinfo.coeff=[];
sol.solverinfo.tmesh=[];

