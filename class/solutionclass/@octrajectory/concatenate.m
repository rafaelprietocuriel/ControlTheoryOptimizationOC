function ocTrj=concatenate(varargin)
ocTrj=[];
if nargin==2
    ocTrj1=varargin{1};
    ocTrj2=varargin{2};

    ocTrj.arcarg=[ocTrj1.arcarg ocTrj2.arcarg];
    ocTrj.x=[ocTrj1.x ocTrj1.x(end)+ocTrj2.x];
    ocTrj.arcinterval=[ocTrj1.arcinterval ocTrj1.arcinterval(end)+ocTrj2.arcinterval(2:end)];
    ocTrj.y=[ocTrj1.y ocTrj2.y];
    ocTrj.timehorizon=ocTrj1.timehorizon+ocTrj2.timehorizon;
    idx=find(diff(ocTrj.x)==0);
    ocTrj.arcposition=[1 idx+1;idx length(ocTrj.x)];

    %[ocTrj.x,ocTrj.arcinterval]=normalizetime(t,ocTrj.arcposition);


    ocTrj=octrajectory(ocTrj);
    didx=find(diff(ocTrj.arcarg)==0);
    if~isempty(didx)
        ocTrj=mergearc(ocTrj,didx);
    end
else
    ocTrj=concatenate(varargin{1},concatenate(varargin{2:end}));
end