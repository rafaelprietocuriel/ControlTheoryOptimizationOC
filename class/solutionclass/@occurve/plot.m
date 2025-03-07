function varargout=plot(varargin)
%
% PLOT basic plot command of an occurve.
%
% PLOT(OCTRJ,XCOORD,YCOORD) plots the coordinates of the occurve
% OCTRJ, where the x-value is given by the 'XCOORD' and the y-value is given by the 'YCOORD'.
% Every coordinate of the occurve is plotted as a separate line. If
% 'XCOORD' or 'YCOORD' is set to 0 the time argument is plotted on the
% corresponding axis.
%
% PLOT(OCTRJ,XCOORD,YCOORD,'PropertyName','PropertyValue', ...)
% 'continuous' can be set to 'on' or 'off'. If it is set to 'on' each of
% the arcs is plotted separately, otherwise not.
%
% PLOT(AHANDLE,OCPT,XCOORD,YCOORD,'PropertyName','PropertyValue', ...)
% plots to the axes with AHANDLE.
%


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

if ~isempty(varargin)
    xcoord=varargin{1};
    if ~isnumeric(xcoord) || any(rem(xcoord,1)) || any(xcoord<0)
        ocmaterror('X coordinate must be a vector of nonnegative integers.')
    end
    varargin(1)=[];
end

if ~isempty(varargin)
    ycoord=varargin{1};
    if isnumeric(ycoord)
        if  any(rem(ycoord,1)) || any(ycoord<0)
            ocmaterror('Y coordinate must be a vector of nonnegative integers.')
        else
            varargin(1)=[];
        end
    else
        xdata='independentvar';
        ycoord=xcoord;
        xcoord=1;
    end
else
    xdata='independentvar';
    ycoord=xcoord;
    xcoord=1;
end
if length(xcoord)>1 && length(ycoord)>1
    ocmaterror('More than one x and y coordinate.')
end
xdataidx=find(strcmpi(varargin,'xdata'));
ydataidx=find(strcmpi(varargin,'ydata'));
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
varargin(remidx)=[];

%newplot
if ~isempty(axhdl)
    h=line(ocCuv,'xcoord',xcoord,'ycoord',ycoord,'xdata',xdata,'ydata',ydata,'Parent',axhdl,varargin{:});
else
    h=line(ocCuv,'xcoord',xcoord,'ycoord',ycoord,'xdata',xdata,'ydata',ydata,varargin{:});
end
if nargout==1
    varargout{1}=h;
end