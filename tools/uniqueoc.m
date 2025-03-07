function [udynPrim,varargout]=uniqueoc(dynPrim,varargin)
%
% UNIQUEOC find unique objects of cell array of oc elements.
%
% UNIQUEOC(DYNPRIM) finds the distinct elements of a cell array of
% dynprimitive objects. Empty cells are removed.
%
% UNIQUEOC(OCELMNTS,OPT) OPT provides options defining relative and
% absolute tolerances describing the same object.
%
% OPT is a oc option structure defining
%       AbsTol  ... absolute tolerance for distinct points (1e-4)
%       RelTol  ... relative tolerance  -"-                (1e-4)
%
% OCPTN = UNIQUEOC(OCELMNTS,...)
%
% [OCPTN,IDX] = UNIQUEOC(OCELMNTS,...), where OCELMNTS(IDX) = OCPTN.

% initialize variables
opt=[];
idx=[];
udynPrim=dynPrim;

if isempty(dynPrim)
    if nargout>=1
        varargout{1}=idx;
    end
    return
end

if isnumeric(dynPrim)
    [udynPrim,idx]=uniquocnumeric(dynPrim,varargin{:});
    if nargout>=2
        varargout{1}=idx;
    end
    return
end
idx=1;
try
    if ~all(cellfun(@isequilibrium,dynPrim))
        ocmatmsg('Input argument is not a cell of Equilibria.')
        return
    end
end
cmpidx=2:length(dynPrim);

if isempty(cmpidx)
    if nargout>=1
        varargout{1}=idx;
    end
    return
end

if nargin>=2
    opt=varargin{1};
end

if isempty(opt)
    opt=defaultocoptions;
end
abstol=getocoptions(opt,'GENERAL','AdmissibleTolerance');
reltol=getocoptions(opt,'GENERAL','AdmissibleTolerance');

while ~isempty(cmpidx)
    remidx=[];
    for ii=1:length(cmpidx)
        dist=reltol*norm(dynPrim{idx(end)})+abstol;

        diffval=distance(dynPrim{idx(end)},dynPrim{cmpidx(ii)});
        if diffval<=dist
            remidx= [remidx ii];
        end
    end
    if isempty(remidx) || remidx(1)~=1
        idx=[idx cmpidx(1)];
    end
    cmpidx([1 remidx])=[];
end

udynPrim=dynPrim(idx);
if nargout>=2
    varargout{1}=idx;
end

function [X,idx]=uniquocnumeric(X,varargin)
opt=[];
idx=1;
cmpidx=2:size(X,2);

if isempty(cmpidx)
    if nargout>=1
        varargout{1}=idx;
    end
    return
end

if nargin>=2
    opt=varargin{1};
end

if isempty(opt)
    opt=defaultocoptions;
end
abstol=getocoptions(opt,'GENERAL','AdmissibleTolerance');
reltol=getocoptions(opt,'GENERAL','AdmissibleTolerance');

while ~isempty(cmpidx)
    remidx=[];
    for ii=1:length(cmpidx)
        dist=reltol*norm(X(:,idx(end)))+abstol;

        diffval=norm(X(:,idx(end))-X(:,cmpidx(ii)));
        if diffval<=dist
            remidx= [remidx ii];
        end
    end
    if isempty(remidx) || remidx(1)~=1
        idx=[idx cmpidx(1)];
    end
    cmpidx([1 remidx])=[];
end

X=X(:,idx);

