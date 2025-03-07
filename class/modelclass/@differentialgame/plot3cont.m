function varargout=plot3cont(varargin)
%
% PLOT3CONT plots the result(s) of a continuation process
%
% PLOT3CONT(OCOBJ,'XDATA',XCOORD,'YDATA',YCOORD,'ZDATA',ZCOORD) plots for
% each continuation process of class 'extremal' stored in OCOBJ the last
% detected solution (entry of the ExtremalSolution field). The values of
% the x/y/z-axis are determined by the coordinate XCOORD/YCOORD/ZCOORD of
% the term specified in 'XDATA'/'YDATA'/'ZDATA'. Possible values for
% 'XDATA', 'YDATA' and 'ZDATA' are:  
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
% The values of XCOORD, YCOORD, ZCOORD specify the coordinates of the
% corresponding terms.
%
% PLOT3CONT(...,'PROPERTYNAME',PROPERTYVALUE,...) 
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
%   'view3d'    : 'on' sets the default three-dimensional view 
%                 'off' uses the view specified by the axis properties.  
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
% PLOT3CONT(AHANDLE,...) plots into the axes with the handle AHANDLE
% instead of into the current axes (gca) 
%
% H=PLOT3CONT(...) returns a column vector of handles to lineseries
% graphics objects, one handle per line. 

axes_hdl=[];
idx=[];
contindex=[];
contfield='';
contclass='';
holdinter=[];
view3Dindex=[];
if ishandle(varargin{1})
    axes_hdl=varargin{1};
    varargin(1)=[];
end
msgstr=nargchk(7,inf,numel(varargin));
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
zvar=varargin{6};
zcoord=varargin{7};

varargin(1:7)=[];

if isempty(ocObj) || isempty(resultstruct) || ~isfield(resultstruct,'Continuation')
    if nargout==1
        varargout{1}=[];
    end
    return
end
contresult=resultstruct.Continuation;

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
        case 'contindex'
            contindex=varargin{ii+1};
            removeidx=[removeidx ii ii+1];
        case 'hold'
            holdinter=strcmpi(varargin{ii+1},'on');
            removeidx=[removeidx ii ii+1];
        case 'view3d'
            view3Dindex=strcmpi(varargin{ii+1},'on');
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
numresult=numel(contresult);
if isempty(idx)
    idx=1:numresult;
end

if min(idx)<1 || max(idx)>numresult
    ocmaterror('Index exceeds number of contiuation results.')
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
if isempty(view3Dindex)
    view3Dindex=0;
end

counter=0;
h=[];
if ~isempty(axes_hdl)
    nextplotori=get(axes_hdl,'NextPlot');
else
    nextplotori=get(gca,'NextPlot');
end
for ii=idx
    if strmatch(contclass,contresult{ii}.ContinuationClassification)
        tagname=[contfield '_' contresult{ii}.ContinuationClassification];
        if ~isfield(contresult{ii},contfield)
            ocmaterror('%s is not a field of the continuation result.')
        end
        counter=counter+1;
        solObj=contresult{ii}.(contfield);
        if isempty(contindex) || numel(idx)>1
            contindex=1:numel(solObj);
        end
        if isstruct(solObj)
            for jj=contindex
                if strcmp(contfield,'ContinuationInformation')
                    switch contclass
                        case {'indifferencesolution','heteroclinic','indifferencedistribution'}
                            currentsolObj=sol2multipath(solObj(jj).data.sol);
                        case 'extremal'
                            currentsolObj=octrajectory(solObj(jj).data.sol);
                    end
                else
                    switch contclass
                        case {'indifferencesolution','heteroclinic','indifferencedistribution'}
                            currentsolObj=sol2multipath(solObj(jj));
                        case {'extremal','limitextremal'}
                            currentsolObj=octrajectory(solObj(jj));
                    end
                end
                %htmp=line(ocTrj,'xdata',xvar,'xcoord',xcoord,'ydata',yvar,'ycoord',ycoord,'zdata',zvar,'zcoord',zcoord,'ocmodel',ocObj,'Tag',[tagname '_Nr:' num2str(ii) '_' num2str(jj)],varargin{:});
                htmp=plot3(currentsolObj,xcoord,ycoord,zcoord,'xdata',xvar,'ydata',yvar,'zdata',zvar,'ocmodel',ocObj,varargin{:});
                set(htmp,'Tag',[tagname '_Nr:' num2str(ii) '_' num2str(jj)]);
                if view3Dindex
                    view(3)
                end
                if ~holdinter && (jj<contindex(end) || ii~=idx(end))
                    figure(gcf)
                    drawnow
                    delete(htmp)
                else
                    if jj==contindex(1)
                        set(gca,'NextPlot','add')
                    end
                    h=[h;htmp];
                    figure(gcf)
                    drawnow
                end
            end
        else
            %htmp=line(solObj,'xdata',xvar,'xcoord',xcoord,'ydata',yvar,'ycoord',ycoord,'zdata',zvar,'zcoord',zcoord,'ocmodel',ocObj,'Tag',[tagname '_Nr:' num2str(ii)],varargin{:});
            htmp=plot3(solObj,xcoord,ycoord,zcoord,'xdata',xvar,'ydata',yvar,'zdata',zvar,'ocmodel',ocObj,varargin{:});
            set(htmp,'Tag',[tagname '_Nr:' num2str(ii)]);
            h=[h;htmp];
        end
    else
        ocmatmsg('Nothing to be plotted!\n')
    end
    if ii==idx(1)
        set(gca,'NextPlot','add')
    end
end
set(gca,'NextPlot',nextplotori)
if view3Dindex
    view(3)
end
if nargout==1
    varargout{1}=h;
end