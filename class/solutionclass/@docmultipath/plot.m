function varargout=plot(varargin)
%
% PLOT basic plot command of a docmultipath.
%

plotarg{1}=[];
if ishandle(varargin{1})
    if ~strcmp(get(varargin{1},'type'),'axes')
        return
    end
    plotarg{1}=varargin{1};
    varargin(1)=[];
end

if ~isdocmultipath(varargin{1})
    ocmaterror('Invalid input argument.')
end
m=multiplicity(varargin{1});
actcolor=[0 0 1;1 0 0;0 0.5 0];
nextplotarg=[];
if ~isempty(plotarg{1}) 
    axh=plotarg{1};
    nextplotarg=get(plotarg{1},'NextPlot');
elseif ~isempty(findobj('Type','axes'))
    nextplotarg=get(gca,'NextPlot');
end
if isempty(plotarg{1}) 
    axh=gca;
end
h=[];
for ii=1:m
    htmp=plot(axh,varargin{1}.solutionclass{ii},varargin{2:end});
    if ~any(strcmp(varargin(2:end),'Color'))
        set(htmp,'Color',actcolor(ii,:))
    end
    if ii==1
        hold on
    end
    h=[h;htmp];
end

if isempty(nextplotarg)
    hold off
else
    set(gca,'NextPlot',nextplotarg);
end
if nargout==1
    varargout{1}=h;
end
