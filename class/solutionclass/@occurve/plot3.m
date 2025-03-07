function varargout=plot3(varargin)
%
% PLOT3 3D plot command for an occurve.
%
% PLOT3(OCTRJ,XCOORD,YCOORD) plots the coordinates XCOORD and YCOORD of
% the dependent variables of the occurve OCTRJ against the independent
% variable. XCOORD and YCOORD have to be single valued..
%
% PLOT3(OCTRJ,XCOORD,YCOORD,'XDATA',XDATA,'YDATA',YDATA,'OCMODEL',OCOBJ) 
%   'XDATA'/'YDATA'  : specifies the plotted terms. Default value is
%                      'dependentvar'. For a detailed list see the LINE
%                      command. 
%   'OCMODEL': OCOBJ is an object of the used model class. These only have
%              to be provided if model specific terms, e.g. state, costate,
%              control, ... are used for 'XDATA' or 'YDATA'.
%  
% PLOT3(AHANDLE,...) plots into the axes with the handle AHANDLE instead of
% into the current axes (gca)
% plots to the axes with AHANDLE.
%  
% PLOT3(...,'CONNECT',CONNECTVAL,...) CONNECTVAL : 0 (default)/1. For 0
% different arcs of the occurve OCTRJ are plotted as separate lines.
% For 1 different arcs are combined into a single line.
% 
%  
% PLOT3(...'PROPERTYNAME',PROPERTYVALUE,...) sets properties to the
% specified property values for all lineseries graphics objects created by
% PLOT3. 
%
% H=PLOT3(...) returns a column vector of handles to lineseries graphics
% objects, one handle per line. 

zcoord=[];
xdata='';
ydata='';
axhdl=[];

if ishandle(varargin{1})
    if ~strcmp(get(varargin{1},'type'),'axes')
        return
    end
    axhdl=varargin{1};
    varargin(1)=[];
end

ocCuv=varargin{1};
varargin(1)=[];
if ~isoccurve(ocCuv)
    ocmaterror('First or second argument is not an ''occurve''.')
end

if isempty(ocCuv)
    if nargout==1
        varargout{1}=[];
    end
    return
end

if length(varargin)<2
    ocmaterror('Not enough input arguments')
end
xcoord=varargin{1};
if ~isnumeric(xcoord) || any(rem(xcoord,1)) || any(xcoord<0)
    ocmaterror('X coordinate must be a vector of nonnegative integers.')
end
varargin(1)=[];

ycoord=varargin{1};
if ~isnumeric(ycoord) || any(rem(ycoord,1)) || any(ycoord<0)
    ocmaterror('Y coordinate must be a vector of nonnegative integers.')
end
varargin(1)=[];

if ~isempty(varargin)
    if isnumeric(varargin{1})
        zcoord=varargin{1};
        if  any(rem(zcoord,1)) || any(zcoord<0)
            ocmaterror('Z coordinate must be a vector of nonnegative integers.')
        else
            varargin(1)=[];
        end
    end
end

xdataidx=find(strcmpi(varargin,'xdata'));
ydataidx=find(strcmpi(varargin,'ydata'));
zdataidx=find(strcmpi(varargin,'zdata'));

colorflag=any(strcmpi(varargin,'color'));
remidx=[];
if isempty(xcoord) && ~isempty(ycoord)
    ocmaterror('Vectors must be the same lengths.')
elseif isempty(xcoord) && isempty(y)
    [varargout{1:nargout}]=plot([],[]);
    return
end

if ~isempty(xdata) && ~isempty(xdataidx)
    ocmatmsg('XData is set to independent variable.')
    varargin(xdataidx:xdataidx+1)=[];
    xdataidx=[];
end
if isempty(xdata) && isempty(xdataidx)
    xdata='dependentvar';
elseif ~isempty(xdataidx)
    xdata=varargin{xdataidx+1};
    remidx=[remidx xdataidx:xdataidx+1];
end
if isempty(ydataidx)
    ydata='dependentvar';
else
    ydata=varargin{ydataidx+1};
    remidx=[remidx ydataidx:ydataidx+1];
end
if isempty(zdataidx)
    zdata='dependentvar';
else
    zdata=varargin{zdataidx+1};
    remidx=[remidx zdataidx:zdataidx+1];
end
varargin(remidx)=[];

if isempty(zcoord)
    zdata=ydata;
    ydata=xdata;
    xdata='independentvar';
    zcoord=ycoord;
    ycoord=xcoord;
    xcoord=1;
end

newplot
if ~isempty(axhdl)
    h=line(ocCuv,'xcoord',xcoord,'ycoord',ycoord,'zcoord',zcoord,'xdata',xdata,'ydata',ydata,'zdata',zdata,'Parent',axhdl,varargin{:});
else
    h=line(ocCuv,'xcoord',xcoord,'ycoord',ycoord,'zcoord',zcoord,'xdata',xdata,'ydata',ydata,'zdata',zdata,varargin{:});

end
if nargout==1
    varargout{1}=h;
end