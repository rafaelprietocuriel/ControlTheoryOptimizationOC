function varargout=pdemesh(ppdeTrj)

hh=pdemesh(points(ppdeTrj),edges(ppdeTrj),triangles(ppdeTrj),ppdeTrj.ppdetrajectory.y(1:end/2));

if nargout==1
    varargout{1}=hh;
end
