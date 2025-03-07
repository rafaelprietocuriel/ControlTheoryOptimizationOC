function varargout=plot3result(varargin)
%
% PLOTCONT plots the result(s) of a continuation process
%
% PLOTCONT(ODEOBJ,'XDATA',XCOORD,'YDATA',YCOORD) plots for each continuation
% process of class 'extremal' stored in ODEOBJ the last detected solution
% (entry of the ExtremalSolution field). The values of the x/y-axis are
% determined by the coordinate XCOORD/YCOORD of the term specified in
% 'XDATA'/'YDATA'. Possible values for 'XDATA' and 'YDATA' are: 
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
% PLOTCONT(...,'PROPERTYNAME',PROPERTYVALUE,...) 
% OCMAT specific pairs 'PROPERTYNAME',PROPERTYVALUE are:
%   'contclass' : 'extremal' continuation of extremal solutions
%                 'indifferencesolution' continuation of indifference
%                                        thresholds 
%   'contfield' : 'ExtremalSolution' last computed solution of a
%                 continuation process  
%                 'ContinuationSolution' solutions computed during the
%                 continuation process.  
%   'index'     : number(s) specifying the continuation process for which
%                 the solution(s) are plotted.
%   'connect'   : (0)/1 see stdocmodel/plot.
% The following values for 'PROPERTYNAME' have only an effect if the
% 'contfield' is set to 'ContinuationSolution'
%   'contindex' : a vector that specifies the indices of solutions that are
%                 plotted.
%   'hold'      : if set to 'on' the axis property 'NextPlot' is set to
%                 'add' during the plotting of the results of the
%                 continuation process. 
%               : if set to 'off' the axis property 'NextPlot' is set to
%                 'replace' during the plotting of the results of the
%                 continuation process. 
%
% For 'PROPERTYNAME' and PROPERTYVALUE usual axes properties can be
% provided.
%
% PLOTCONT(AHANDLE,...) plots into the axes with the handle AHANDLE instead
% of into the current axes (gca)
%
% H=PLOTCONT(...) returns a column vector of handles to lineseries graphics
% objects, one handle per line. 
axes_hdl=[];
idx=[];
contfield='';
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
resultStruct=result(odeObj);
xvar=varargin{2};
xcoord=varargin{3};
yvar=varargin{4};
ycoord=varargin{5};
zvar=varargin{6};
zcoord=varargin{7};

varargin(1:7)=[];

if isempty(resultStruct)
    if nargout==1
        varargout{1}=[];
    end
    return
end
if mod(numel(varargin),2)
    ocmaterror('Input arguments have to be called in name/property value pairs.')
end

removeidx=[];
for ii=1:2:numel(varargin)
    switch lower(varargin{ii})
        case 'index'
            idx=varargin{ii+1};
            removeidx=[removeidx ii ii+1];
        case 'contfield'
            contfield=varargin{ii+1};
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
if isempty(contfield)
    contfield='octrajectory';
end

counter=0;
h=[];
if ~isempty(axes_hdl)
    nextplotori=get(axes_hdl,'NextPlot');
else
    nextplotori=get(gca,'NextPlot');
end

if ~isfield(resultStruct,contfield)
    return
end
numresult=numel(resultStruct.(contfield));

if isempty(idx)
    idx=1:numresult;
end

if min(idx)<1 || max(idx)>numresult
    ocmaterror('Index exceeds number of continuation results.')
end
tagname=contfield;
for ii=idx
    counter=counter+1;
    htmp=plot3(resultStruct.(contfield){ii},xcoord,ycoord,zcoord,'xdata',xvar,'ydata',yvar,'zdata',zvar,'ocmodel',odeObj,varargin{:});
    set(htmp,'Tag',[tagname '_Nr:' num2str(ii)])
    h=[h;htmp];
    if ii==idx(1)
        set(gca,'NextPlot','add')
    end
end
set(gca,'NextPlot',nextplotori)
if nargout==1
    varargout{1}=h;
end
