function ppdeTrjNew=spacedeval(ppdeTrj,znew,ppdeObj,varargin)

ppdeTrjNew=ppdeTrj;
y=ppdeTrj.octrajectory.y;
femdat=femdata(ppdeTrj);
z0=femdat.grid.x;
scoord=spatialcoordinates(ppdeObj,ppdeTrj);
nonscoord=1:size(y,1);
nonscoord(scoord(:))=[];

if length(z0)>=length(znew)
    nanval=repmat(NaN,length(z0)-length(znew),size(y,2));
else
end
ynew=zeros(size(scoord,2)*length(znew)+length(nonscoord),size(y,2));
for ii=1:size(scoord,2)
    ytmp=interp1(z0,y(scoord(:,ii),:),znew,'linear');
    if size(y,2)==1
        ytmp=ytmp(:);
    end
    ynew(scoord(:,ii),:)=[ytmp;nanval];
end
ynew(nonscoord,:)=y(nonscoord,:);
ynew(isnan(ynew(:,1)),:)=[];

ppdeTrjNew.octrajectory.y=ynew;
femdatNew=femdat;
femdatNew.grid=Interval(znew);
femdatNew.gridnum=length(znew);
ppdeTrjNew.femdata=femdatNew;