function varargout=plotlimitset(varargin)
%
% PLOTLIMITSET plots the stored limitsets
%
% PLOTLIMITSET(OCOBJ,'XDATA',XCOORD,'YDATA',YCOORD) plots the equilibria
% stored in OCOBJ. The values of the x/y-axis are determined by the
% coordinate XCOORD/YCOORD of the term specified in 'XDATA'/'YDATA'.
% Possible values for 'XDATA' and 'YDATA' are:  
%
%   'independentvar'
%   'time'
%   'dependentvar'
%   'state'
%   'costate',
%   'control'
%   'lagrangemultiplier'
%   'canonicalsystem'
%   'hamiltonian'
%   'userfunction'
% 
% The values of XCOORD and YCOORD specify the coordinates of the
% corresponding terms. Either XCOORD or YCOORD can be a vector.
%
% PLOTLIMITSET(...,'PROPERTYNAME',PROPERTYVALUE,...) 
% OCMAT specific pairs 'PROPERTYNAME',PROPERTYVALUE are:
%   'limitclass': 'equilibrium' the equilibria stored in the result field
%                 'Equilibrium'
%   'index'     : integer vector specifying the index of the limitsets that
%                 are plotted. 
%
% For 'PROPERTYNAME' and PROPERTYVALUE usual axes properties can be
% provided.
%
% PLOTLIMITSET(AHANDLE,...) plots into the axes with the handle AHANDLE
% instead of into the current axes (gca)
%
% H=PLOTLIMITSET(...) returns a column vector of handles to lineseries
% graphics objects, one handle per line. 

axes_hdl=[];
idx=[];
limitclass='';
showtag=0;
if ishandle(varargin{1})
    axes_hdl=varargin{1};
    varargin(1)=[];
end
msgstr=nargchk(5,inf,numel(varargin));
if ~isempty(msgstr)
    ocmaterror(msgstr);
end
msgstr=nargoutchk(0,1,nargout);
if ~isempty(msgstr)
    ocmaterror(msgstr);
end

if ~isocsolutionclass(varargin{1})
    ocmaterror('Inputargument %s is not an ocmodel.',inputname(1));
end

odeObj=varargin{1};
xvar=varargin{2};
xcoord=varargin{3};
yvar=varargin{4};
ycoord=varargin{5};

varargin(1:5)=[];
if mod(numel(varargin),2)
    ocmaterror('Input arguments have to be called in name/property value pairs.')
end

removeidx=[];
for ii=1:2:numel(varargin)
    switch lower(varargin{ii})
        case 'index'
            idx=varargin{ii+1};
            removeidx=[removeidx ii ii+1];
        case 'limitclass'
            limitclass=varargin{ii+1};
            removeidx=[removeidx ii ii+1];
        case 'showtag'
            showtag=varargin{ii+1};
            removeidx=[removeidx ii ii+1];
    end
end
if idx==0
    return
end
varargin(removeidx)=[];
if ~isempty(axes_hdl)
    varargin{end+1}='Parent';
    varargin{end+1}=axes_hdl;
end
if isempty(limitclass)
    limitclass='Equilibrium';
end
if isempty(showtag)
    showtag=0;
end

if strcmpi(limitclass,'Equilibrium')
    limitObj=equilibrium(odeObj);
elseif strcmpi(limitclass,'PeriodicSolution')
    limitObj=periodicsolution(odeObj);
end

if isempty(limitObj)
    return
end
numresult=numel(limitObj);
if isempty(idx)
    idx=1:numresult;
end
if min(idx)<1 || max(idx)>numresult
    ocmaterror('Index exceeds number of continuation results.')
end

counter=0;
h=[];
tagname=limitclass;
if ~isempty(axes_hdl)
    nextplotori=get(axes_hdl,'NextPlot');
else
    nextplotori=get(gca,'NextPlot');
end
for ii=idx
    counter=counter+1;
    solObj=limitObj{ii};
    h(counter)=plot(solObj,xcoord,ycoord,'xdata',xvar,'ydata',yvar,'ocmodel',odeObj,varargin{:});
    set(h(counter),'Tag',[tagname '_Nr:' num2str(ii)])
    if ii==idx(1)
        set(gca,'NextPlot','add')
    end
    if showtag
        xdata=get(h(counter),'XData');
        ydata=get(h(counter),'YData');
        htext=text(xdata(1),ydata(1),num2str(ii));
        set(htext,'Tag',[tagname '_Nr:' num2str(ii) '_' num2str(ii)],'Interpreter','latex')
    end
end
set(gca,'NextPlot',nextplotori)
if nargout==1
    varargout{1}=h;
end
