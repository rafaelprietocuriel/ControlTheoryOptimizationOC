function varargout=plotmatcont(varargin)


axes_hdl=[];
idx=[];
contfield='';
contclass='';
holdinter=[];
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
resultstruct=result(ocObj);
xvar=varargin{2};
xcoord=varargin{3};
yvar=varargin{4};
ycoord=varargin{5};

varargin(1:5)=[];

if isempty(ocObj) || isempty(resultstruct) || ~isfield(resultstruct,'MatContContinuation')
    if nargout==1
        varargout{1}=[];
    end
    return
end
contresult=resultstruct.MatContContinuation;

if mod(numel(varargin),2)
    ocmaterror('Input arguments have to be called in name/property value pairs.')
end

% possible properties/beside native MATLAB plot properties
% Index ... index of continuation elements
% Contfield ... field of 
% Contclass ... 'extremal2ep'
% Holdcontsolution ... hold 
removeidx=[];
for ii=1:2:numel(varargin)
    switch lower(varargin{ii})
        case 'index'
            idx=varargin{ii+1};
            removeidx=[removeidx ii ii+1];
        case 'contfield'
            contfield=varargin{ii+1};
            removeidx=[removeidx ii ii+1];
        case 'contclass'
            contclass=varargin{ii+1};
            removeidx=[removeidx ii ii+1];
        case 'hold'
            holdinter=strcmpi(varargin{ii+1},'on');
            removeidx=[removeidx ii ii+1];
    end
end
varargin(removeidx)=[];
if ~isempty(axes_hdl)
    varargin{end+1}='Parent';
    varargin{end+1}=axes_hdl;
end
numresult=numel(contresult);
if isempty(idx)
    idx=1:numresult;
end

if min(idx)<1 || max(idx)>numresult
    ocmaterror('Index exceeds number of contiuation results.')
end
if isempty(contclass)
    contclass='modelequilibrium';
end
if isempty(contfield)
    contfield='ContinuationSolution';
end
if isempty(holdinter)
    holdinter=0;
end

counter=0;
h=[];
tagname=[contfield '_' contclass];
for ii=idx
    if strcmpi(contresult{ii}.ContinuationClassification,contclass)
        if ~isfield(contresult{ii},contfield)
            ocmaterror('%s is not a field of the continuation result.')
        end
        counter=counter+1;
        solObj=contresult{ii}.(contfield);
        numObj=numel(solObj);
        if isstruct(solObj) && numObj>1
            for jj=1:numObj
                if strcmp(contfield,'ContinuationInformation')
                    ocTrj=octrajectory(solObj(jj).data.sol);
                else
                    ocTrj=octrajectory(solObj(jj));
                end
                htmp=line(ocTrj,'xdata',xvar,'xcoord',xcoord,'ydata',yvar,'ycoord',ycoord,'ocmodel',ocObj,'Tag',[tagname '_Nr:' num2str(ii) '_' num2str(jj)],varargin{:});
%                 figure(gcf)
%                 drawnow
                if ~holdinter && jj<numObj
                    delete(htmp)
                else
                    h=[h htmp];
                end
            end
        else
            if isstruct(solObj)
                solObj=octrajectory(solObj);
            end
            h(counter)=line(solObj,'xdata',xvar,'xcoord',xcoord,'ydata',yvar,'ycoord',ycoord,'ocmodel',ocObj,'Tag',[tagname '_Nr:' num2str(ii)],varargin{:});
        end
    end
end
if nargout==1
    varargout{1}=h;
end
