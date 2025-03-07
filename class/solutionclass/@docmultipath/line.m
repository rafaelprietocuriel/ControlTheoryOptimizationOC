function varargout=line(docMultiPath,varargin)
%
% LINE basic plot command of an octrajectory.
%
% LINE(OCTRJ,'XVAR',XVAR,'XCOORD',XCOORD,'YVAR',YVAR,'YCOORD',YCOORD,'OCMODEL',OCOBJ)
%
% LINE(OCTRJ,XCOORD,YCOORD)
%
% LINE(OCTRJ,XCOORD,YCOORD,ZCOORD)
%
% LINE(OCTRJ,XCOORD,YCOORD,ZCOORD,'PropertyName',propertyvalue,...)
%
% further properties are: 'OCModel'
%           'Connect' ... determines if line is plotted separately for
%           every arc 

h=[];

m=multiplicity(docMultiPath);
for ii=1:m
    htmp=line(docMultiPath.solutionclass{ii},varargin{1:end});
    h=[h;htmp];
end

if nargout==1
    varargout{1}=h;
end