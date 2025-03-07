function varargout=plot(varargin)
%
% PLOT plot command for an octrajectory.
%
% PLOT(DOCTRJ,COORD) plots the coordinates COORD (integer vector) of the
% dependent variables of the octrajectory DOCTRJ against the independent
% variable.
%
% PLOT(DOCTRJ,XCOORD,YCOORD) plots the coordinates XCOORD and YCOORD of
% the dependent variables of the octrajectory DOCTRJ. XCOORD or YCOORD can
% be an integer vector.
%
% PLOT(DOCTRJ,COORD,'YDATA',YDATA,'OCMODEL',OCOBJ) 
%   'YDATA'  : specifies the plotted terms. Default value is
%              'dependentvar'. For a detailed list see the LINE command.
%   'OCMODEL': OCOBJ is an object of the used model class. These only have
%              to be provided if model specific terms, e.g. state, costate,
%              control, ... are used for 'YDATA'
%
% PLOT(DOCTRJ,XCOORD,YCOORD,'XDATA',XDATA,'YDATA',YDATA,'OCMODEL',OCOBJ) 
%   'XDATA'/'YDATA'  : specifies the plotted terms. Default value is
%                      'dependentvar'. For a detailed list see the LINE
%                      command. 
%   'OCMODEL': OCOBJ is an object of the used model class. These only have
%              to be provided if model specific terms, e.g. state, costate,
%              control, ... are used for 'XDATA' or 'YDATA'.
%  
% PLOT(AHANDLE,...) plots into the axes with the handle AHANDLE instead of
% into the current axes (gca)
%  
% PLOT(...,'CONNECT',CONNECTVAL,...) CONNECTVAL : 0 (default)|1. For 0
% different arcs of the octrajectory DOCTRJ are plotted as separate lines.
% For 1 different arcs are combined into a single line.
% 
%  
% PLOT(...'PROPERTYNAME',PROPERTYVALUE,...) sets properties to the
% specified property values for all lineseries graphics objects created by
% PLOT. 
%
% H=PLOT(...) returns a column vector of handles to lineseries graphics
% objects, one handle per line. 

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

docTrj=varargin{1};
varargin(1)=[];
if ~isdoctrajectory(docTrj)
    ocmaterror('First or second argument is not a ''doctrajectory''.')
end

if isempty(docTrj)
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
    h=line(docTrj,'xcoord',xcoord,'ycoord',ycoord,'xdata',xdata,'ydata',ydata,'Parent',axhdl,'LineStyle','none','Marker','x',varargin{:});
else
    h=line(docTrj,'xcoord',xcoord,'ycoord',ycoord,'xdata',xdata,'ydata',ydata,'LineStyle','none','Marker','x',varargin{:});
end
if nargout==1
    varargout{1}=h;
end

if ~colorflag
    colour=get(gca,'ColorOrder');
    numcolour=size(colour,1);
    countercolor=0;
    while ~isempty(h)
        tagflag=get(h(1),'Tag');
        tagnum=str2double(tagflag)+1;
        if isempty(tagflag) && 1<=tagnum && tagnum<=numcolour
            if countercolor<numcolour
                countercolor=countercolor+1;
            else
                countercolor=1;
            end
            set(h(1),'Color',colour(rem(countercolor,numcolour)+1,:))
        else
            if ~isempty(tagflag)
                set(h(1),'Color',colour(rem(tagnum,numcolour)+1,:))
            end
        end
        h(1)=[];
    end
end