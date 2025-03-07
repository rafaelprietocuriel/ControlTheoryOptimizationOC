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
resultname='';
phaseshift={[]};
showspp=0;
showflat=[];
showtag=0;
discrete=[];
initialpoint=[];

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
        case 'phaseshift'
            phaseshift=varargin{ii+1};
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
        case 'resultname'
            resultname=varargin{ii+1};
            removeidx=[removeidx ii ii+1];
        case 'discrete'
            discrete=varargin{ii+1};
            removeidx=[removeidx ii ii+1];
        case 'initialpoint'
            initialpoint=varargin{ii+1};
            removeidx=[removeidx ii ii+1];
    end
end
if nargout==1
    varargout{1}=[];
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
if isempty(resultname)
    resultname='Equilibrium';
end
if isempty(discrete)
    discrete=0;
end

if strcmpi(limitclass,'Equilibrium')
    %limitObj=equilibrium(ocObj);
    limitObj=result(ocObj,resultname);
elseif strcmpi(limitclass,'PeriodicSolution') || strcmpi(limitclass,'LimitCycle')
    limitObj=periodicsolution(ocObj);
    if isnumeric(phaseshift)
        limitObj{1}=shiftphase(limitObj{1},phaseshift);
    end
else
    limitObj=result(ocObj,limitclass);        
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
    htmp=plot(solObj,xcoord,ycoord,'xdata',xvar,'ydata',yvar,'ocmodel',ocObj,varargin{:});
    if initialpoint
        xdata=get(htmp(1),'XData');
        ydata=get(htmp(1),'YData');
        hinit=plot(xdata(1),ydata(1));
        set(hinit,'Color',[1 0 0],'Marker','x','MarkerSize',10)
    end
    if discrete
        htmp2=zeros(length(independentvar(solObj)),1);
        counter=0;
        for jj=1:length(htmp)
            xdata=get(htmp(jj),'XData');
            ydata=get(htmp(jj),'YData');
            clr=get(htmp(jj),'Color');
            nextplot=get(gca,'NextPlot');
            set(gca,'NextPlot','add');
            for kk=1:length(xdata)
                counter=counter+1;
                htmp2(kk)=plot(xdata(kk),ydata(kk),'Marker','x','Color',clr);
                set(htmp2(kk),'Tag',[tagname '_Nr:' num2str(ii) '_' num2str(jj) '_Pt:' num2str(counter)])
            end
        end
        delete(htmp)
        set(gca,'NextPlot',nextplot);
        h=[h(:);htmp2(:)];
    else
        set(htmp,'Tag',[tagname '_Nr:' num2str(ii) '_' num2str(ii)])
        h=[h;htmp];
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
                xdata=get(htmp,'XData');
                ydata=get(htmp,'YData');
                [b diff]=isspp(ocObj,solObj);
                htext=text(xdata(1),ydata(1),sprintf('$%s_{(%d)}$',num2str(ii),diff));
                set(htext,'Tag',[tagname '_Nr:' num2str(ii) '_' num2str(ii)],'Interpreter','latex')
            else
                xdata=get(htmp,'XData');
                ydata=get(htmp,'YData');
                htext=text(xdata(1),ydata(1),num2str(ii));
                set(htext,'Tag',[tagname '_Nr:' num2str(ii) '_' num2str(ii)],'Interpreter','latex')
            end
        end
        if showspp && ~showtag
            xdata=get(htmp,'XData');
            ydata=get(htmp,'YData');
            [b diff]=isspp(ocObj,solObj);
            htext=text(xdata(1),ydata(1),sprintf('$%s_{(%d)}$',form,diff));
            set(htext,'Tag',[tagname '_Nr:' num2str(ii) '_' num2str(ii)],'Interpreter','latex')
        elseif showflat
            xdata=get(htmp,'XData');
            ydata=get(htmp,'YData');
            htext=text(xdata(1),ydata(1),form);
            set(htext,'Tag',[tagname '_Nr:' num2str(ii) '_' num2str(ii)])
        end
    end
end
set(gca,'NextPlot',nextplotori)
if nargout==1
    varargout{1}=h;
end
