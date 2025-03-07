function varargout=plotcont(varargin)
%
% PLOTCONT plots the result(s) of a continuation process
%
% PLOTCONT(OCOBJ,'XDATA',XCOORD,'YDATA',YCOORD) plots for each continuation
% process of class 'extremal' stored in OCOBJ the last detected solution
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
contindex=[];
showtag=[];
contfield='';
contclass='';
contname='';
discrete=[];
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
        case 'contfield'
            contfield=varargin{ii+1};
            removeidx=[removeidx ii ii+1];
        case 'contclass'
            contclass=varargin{ii+1};
            removeidx=[removeidx ii ii+1];
        case 'contindex'
            contindex=varargin{ii+1};
            removeidx=[removeidx ii ii+1];
        case 'contname'
            contname=varargin{ii+1};
            removeidx=[removeidx ii ii+1];
        case 'hold'
            holdinter=strcmpi(varargin{ii+1},'on');
            removeidx=[removeidx ii ii+1];
        case 'showtag'
            showtag=varargin{ii+1};
            removeidx=[removeidx ii ii+1];
        case 'discrete'
            discrete=varargin{ii+1};
            removeidx=[removeidx ii ii+1];
    end
end
if idx==0
    if nargout==1
        varargout{1}=[];
    end
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
if isempty(discrete)
    discrete=0;
end
if isempty(contname)
    contresultStruct=contresult(ocObj);
else
    contresultStruct=result(ocObj,contname);
end
if isempty(contresultStruct)
    if nargout==1
        varargout{1}=[];
    end
    return
end
numresult=numel(contresultStruct);
if isempty(idx)
    idx=1:numresult;
end

if min(idx)<1 || max(idx)>numresult
    ocmaterror('Index exceeds number of continuation results.')
end
if isempty(contclass)
    contclass='extremal';
end
if isempty(contfield)
    contfield='ExtremalSolution';
end
if isempty(holdinter)
    holdinter=0;
end

counter=0;
h=[];
if ~isempty(axes_hdl)
    nextplotori=get(axes_hdl,'NextPlot');
else
    nextplotori=get(gca,'NextPlot');
end
for ii=idx
    if strmatch(contclass,contresultStruct{ii}.ContinuationClassification) & isfield(contresultStruct{ii},contfield)
        tagname=[contfield '_' contresultStruct{ii}.ContinuationClassification];
%         if ~isfield(contresultStruct{ii},contfield)
%             ocmaterror('%s is not a field of the continuation result.')
%         end
        counter=counter+1;
        solObj=contresultStruct{ii}.(contfield);
        if isempty(contindex) || numel(idx)>1
            contindex=1:numel(solObj);
        end
        if isstruct(solObj)
            for jj=contindex
                if strcmp(contfield,'ContinuationInformation')
                    switch contclass
                        case {'extremal'}
                            currentsolObj=sol2mmultipath(solObj(jj).data.sol,ocObj);
                    end
                else
                    switch contclass
                        case {'extremal'}
                            currentsolObj=sol2mmultipath(solObj(jj),ocObj);
                        case {'indifferencesolution'}
                            currentsolObj=mmultipath(solObj(jj),ocObj);
                    end
                end
                for kk=1:numberofparts(currentsolObj)
                    htmp=plot(currentsolObj(kk),xcoord,ycoord,'xdata',xvar,'ydata',yvar,'ocmodel',ocObj.Model{kk},varargin{:});
                    if length(htmp)==1
                        set(htmp,'Tag',[contname '_' tagname '_Part:' num2str(kk) '_Nr:' num2str(ii) '_' num2str(jj)])
                    else
                        for ll=1:length(htmp)
                            set(htmp(ll),'Tag',[contname '_' tagname '_Part:' num2str(kk) '_Nr:' num2str(ii) '_' num2str(jj) '_' num2str(ll)])
                        end
                    end
                    if holdinter
                        hold on
                        h=[h;htmp];
                    elseif jj~=contindex(end)
                        delete(htmp);
                    else
                        h=[h;htmp];
                    end
                end
                figure(gcf)
                %pause(0.1)
            end
        else
            if iscell(solObj) || ismmultipath(solObj)
                if iscell(solObj)
                    solObj=mmultipath(solObj);
                    par=parametervalue(ocObj);
                    ocObj=multimodel([submodelname(ocObj) submodelname(ocObj)],modelname(ocObj));
                    ocObj=changeparametervalue(ocObj,[par,par]);
                end
                for jj=1:numberofparts(solObj)
                    htmp=plot(solObj(jj),xcoord,ycoord,'xdata',xvar,'ydata',yvar,'ocmodel',ocObj.Model{jj},varargin{:});
                    if  discrete
                        htmp2=zeros(length(independentvar(solObj(jj))),1);
                        counter=0;
                        for kk=1:length(htmp)
                            xdata=get(htmp(kk),'XData');
                            ydata=get(htmp(kk),'YData');
                            clr=get(htmp(kk),'Color');
                            nextplot=get(gca,'NextPlot');
                            set(gca,'NextPlot','add');
                            for ll=1:length(xdata)
                                counter=counter+1;
                                htmp2(ll)=plot(xdata(ll),ydata(ll),'Marker','x','Color',clr);
                                set(htmp2(ll),'Tag',[contname '_' tagname '_Nr:' num2str(ii) '_' num2str(ll) '_Pt:' num2str(counter)])
                            end
                        end
                        delete(htmp)
                        set(gca,'NextPlot',nextplot);
                        h=[h(:);htmp2(:)];
                    else
                        if length(htmp)==1
                            set(htmp,'Tag',[contname '_' tagname '_Part:' num2str(jj) '_Nr:' num2str(ii)])
                        else
                            for kk=1:length(htmp)
                                set(htmp(kk),'Tag',[contname '_' tagname '_Part:' num2str(jj) '_Nr:' num2str(ii) '_' num2str(kk)])
                            end
                        end
                        h=[h;htmp];
                    end
                end
            end
        end
    end
    if ii==idx(1)
        set(gca,'NextPlot','add')
    end
end
set(gca,'NextPlot',nextplotori)
if nargout==1
    varargout{1}=h;
end
