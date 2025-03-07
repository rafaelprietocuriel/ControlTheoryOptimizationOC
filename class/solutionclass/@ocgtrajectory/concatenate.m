function ocTrj=concatenate(varargin)
ocTrj=[];
if nargin==2
    ocTrj1=ocgtrajectory(varargin{1});
    ocTrj2=ocgtrajectory(varargin{2});

    ocTrj.arcarg=[ocTrj1.octrajectory.arcarg ocTrj2.octrajectory.arcarg];
    ocTrj.x=[ocTrj1.octrajectory.x ocTrj1.octrajectory.x(end)+ocTrj2.octrajectory.x];
    ocTrj.arcinterval=[ocTrj1.octrajectory.arcinterval ocTrj1.octrajectory.arcinterval(end)+ocTrj2.octrajectory.arcinterval(2:end)];
    ocTrj.y=[ocTrj1.octrajectory.y ocTrj2.octrajectory.y];
    ocTrj.timehorizon=ocTrj1.octrajectory.timehorizon+ocTrj2.octrajectory.timehorizon;
    idx=find(diff(ocTrj.x)==0);
    ocTrj.arcposition=[1 idx+1;idx length(ocTrj.x)];
    ocTrj.odenum=[ocTrj1.odenum ocTrj1.odenum];
    %[ocTrj.x,ocTrj.arcinterval]=normalizetime(t,ocTrj.arcposition);


    didx=find(diff(ocTrj.arcarg)==0);
    ocTrj=ocgtrajectory(ocTrj);
    if~isempty(didx)
        ocTrj=mergearc(ocTrj,didx);
    end
else
    ocTrj=concatenate(varargin{1},concatenate(varargin{2:end}));
end