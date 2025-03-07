function varargout=plot3limitset(varargin)
%
% PLOT3LIMITSET plots the stored limitsets
%
% PLOT3LIMITSET(OCOBJ,'XDATA',XCOORD,'YDATA',YCOORD) plots the equilibria
% stored in OCOBJ. The values of the x/y/z-axis are determined by the
% coordinate XCOORD/YCOORD/ZCOORD of the term specified in
% 'XDATA'/'YDATA'/'ZDATA'. Possible values for 'XDATA', 'YDATA' and 'ZDATA' 
% are:   
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
% The values of XCOORD, YCOORD and ZCOORD specify the coordinates of the
% corresponding terms.
%
% PLOT3LIMITSET(...,'PROPERTYNAME',PROPERTYVALUE,...) 
% OCMAT specific pairs 'PROPERTYNAME',PROPERTYVALUE are:
%   'limitclass': 'equilibrium' the equilibria stored in the result field
%                 'Equilibrium'
%   'index'     : integer vector specifying the index of the limitsets that
%                 are plotted. 
%   'view3d'    : 'on' sets the default three-dimensional view 
%                 'off' uses the view specified by the axis properties.  
%
% For 'PROPERTYNAME' and PROPERTYVALUE usual axes properties can be
% provided.
%
% PLOT3LIMITSET(AHANDLE,...) plots into the axes with the handle AHANDLE
% instead of into the current axes (gca)
%
% H=PLOT3LIMITSET(...) returns a column vector of handles to lineseries
% graphics objects, one handle per line. 

axes_hdl=[];
idx=[];
limitclass='';
view3Dindex=[];
showspp=0;
showflat=[];
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

ocObj=varargin{1};
xvar=varargin{2};
xcoord=varargin{3};
yvar=varargin{4};
ycoord=varargin{5};
zvar=varargin{6};
zcoord=varargin{7};

varargin(1:7)=[];
if mod(numel(varargin),2)
    ocmaterror('Input arguments have to be called in name/property value pairs.')
end

% possible properties/beside native MATLAB plot properties
% Index ... index of continuation elements
% Contfield ... field of 
% Contclass ... 'extremal2ep'
removeidx=[];
for ii=1:2:numel(varargin)
    switch lower(varargin{ii})
        case 'index'
            idx=varargin{ii+1};
            removeidx=[removeidx ii ii+1];
        case 'limitclass'
            limitclass=varargin{ii+1};
            removeidx=[removeidx ii ii+1];
        case 'view3d'
            view3Dindex=strcmpi(varargin{ii+1},'on');
            removeidx=[removeidx ii ii+1];
        case 'showspp'
            showspp=varargin{ii+1};
            removeidx=[removeidx ii ii+1];
        case 'showflat'
            showflat=varargin{ii+1};
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
if isempty(showtag)
    showtag=0;
end
if isempty(limitclass)
    limitclass='Equilibrium';
end

if strcmp(limitclass,'Equilibrium')
    limitObj=equilibrium(ocObj);
elseif strcmp(limitclass,'PeriodicSolution') || strcmp(limitclass,'LimitCycle')
    limitObj=periodicsolution(ocObj);
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
    htmp=plot3(solObj,xcoord,ycoord,zcoord,'xdata',xvar,'ydata',yvar,'zdata',zvar,'ocmodel',ocObj,varargin{:});
    set(htmp,'Tag',[tagname '_Nr:' num2str(ii) '_' num2str(ii)])
    h=[h;htmp];
    if view3Dindex
        view(3)
    end
    if ii==idx(1)
        set(gca,'NextPlot','add')
    end
    if showflat
        if isflat(ocObj,solObj)
            form='F';
        else
            form='H';
        end
    else
        form='';
    end
    if showtag
        if showspp
            if isspp(ocObj,solObj)
                xdata=get(htmp,'XData');
                ydata=get(htmp,'YData');
                zdata=get(htmp,'ZData');
                text(xdata(1),ydata(1),zdata(1),'SPP')
            end
        else
            xdata=get(htmp,'XData');
            ydata=get(htmp,'YData');
            zdata=get(htmp,'ZData');
            htext=text(xdata(1),ydata(1),zdata(1),num2str(ii));
            set(htext,'Tag',[tagname '_Nr:' num2str(ii) '_' num2str(ii)],'Interpreter','latex')
        end
    end
    if showspp && ~showtag
        xdata=get(htmp,'XData');
        ydata=get(htmp,'YData');
        zdata=get(htmp,'ZData');
        [b diff]=isspp(ocObj,solObj);
        htext=text(xdata(1),ydata(1),zdata(1),sprintf('$%s_{(%d)}$',form,diff));
        set(htext,'Tag',[tagname '_Nr:' num2str(ii) '_' num2str(ii)],'Interpreter','latex')
    elseif showflat
        xdata=get(htmp,'XData');
        ydata=get(htmp,'YData');
        zdata=get(htmp,'ZData');
        htext=text(xdata(1),ydata(1),zdata(1),form);
        set(htext,'Tag',[tagname '_Nr:' num2str(ii) '_' num2str(ii)])
    end
end
set(gca,'NextPlot',nextplotori)
if nargout==1
    varargout{1}=h;
end
