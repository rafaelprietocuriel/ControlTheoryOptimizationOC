function varargout=pdeplot(ppdeTrj)

hh=pdeplot(points(ppdeTrj),edges(ppdeTrj),triangles(ppdeTrj),'xydata',ppdeTrj.ppdetrajectory.y);

if nargout==1
    varargout{1}=hh;
end
