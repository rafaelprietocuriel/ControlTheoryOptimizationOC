function varargout=plotresult(varargin)
%
% PLOTRESULT plots the result(s) of a continuation process
%


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

varargin(1:5)=[];

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
    htmp=plot(resultStruct.(contfield){ii},xcoord,ycoord,'xdata',xvar,'ydata',yvar,'ocmodel',odeObj,varargin{:});
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
