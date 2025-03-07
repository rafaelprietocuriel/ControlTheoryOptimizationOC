function dynPrim=dynprimitive(varargin)
%
%
idx=find(cellfun('isempty',varargin));
if ~isempty(idx)
    varargin(idx)=[];
    dynPrim=dynprimitive(varargin{:});
    return
end

switch nargin
    case 0
        dynPrim.period=[];
        ocTrj=octrajectory([]);
        dynPrim=class(dynPrim,'dynprimitive',ocTrj);
        
    case 1
        if isstruct(varargin{1})
            [ocTrj,period]=adddefault(varargin{1});
            dynPrim.period=period;
            dynPrim=class(dynPrim,'dynprimitive',octrajectory(ocTrj));
        elseif isdynprimitive(varargin{1})
            dynPrim=varargin{1};
        elseif isocasymptotic(varargin{1})
            dynPrimStruct.octrajectory=struct(octrajectory(varargin{1}));
            dynPrimStruct.octrajectory.solverinfo=[];
            dynPrimStruct.period=dynPrimStruct.octrajectory.arcinterval(end);
            dynPrim=dynprimitive(dynPrimStruct);
        elseif isoctrajectory(varargin{1})
            dynPrimStruct=struct(varargin{1});
            dynPrimStruct.period=dynPrimStruct.arcinterval(end);
            dynPrim=dynprimitive(dynPrimStruct);
        end
    case 2
        idx=find(cellfun('isclass',varargin,'stdocmodel')|cellfun('isclass',varargin,'differentialgame'));
        if ~isempty(idx)
            ocObj=varargin{idx};
        else
            ocObj=[];
        end
        varargin(idx)=[];
        idx=find(cellfun('isclass',varargin,'double'));
        if length(idx)==2
            if numel(varargin{idx(1)})==1
                arcarg=varargin{idx(1)};
                varargin(idx(1))=[];
            elseif numel(varargin{idx(2)}==1)
                arcarg=varargin{idx(2)};
                varargin(idx(2))=[];
            else
                error('')
            end
            y=varargin{1};
            ocTrj.arcarg=arcarg;
            ocTrj.y=y(:);
            varargin=[];
        else
            ocTrj=[];
        end
        if ~isempty(varargin)
            ocTrj=varargin{1};
        end
        if ~isempty(ocObj)
            ocTrj.modelname=modelname(ocObj);
            ocTrj.modelparameter=parametervalue(ocObj);
        end
        dynPrim=dynprimitive(ocTrj);
    case 3 % equilibrium y, arcarg, linearization
        idx=find(cellfun('isclass',varargin,'double'));
        if length(idx)==3
            if numel(varargin{1})==1
                arcarg=varargin{1};
                varargin(1)=[];
            elseif numel(varargin{2}==1)
                arcarg=varargin{2};
                varargin(2)=[];
            elseif numel(varargin{3}==1)
                arcarg=varargin{3};
                varargin(3)=[];
            else
                error('')
            end
            if prod(size(varargin{1}))==numel(varargin{1})
                y=varargin{1};
                varargin(1)=[];
            elseif prod(size(varargin{2}))==numel(varargin{2})
                y=varargin{2};
                varargin(2)=[];
            else
                error('')
            end
            J=varargin{1};
            ocTrj.arcarg=arcarg;
            ocTrj.y=y(:);
            ocTrj.linearization=J;
        else
            error('')
        end
        dynPrim=dynprimitive(ocTrj);
    case 4 % ocObj, equilibrium y, arcarg, linearization
        idx=find(cellfun('isclass',varargin,'stdocmodel')|cellfun('isclass',varargin,'differentialgame'));
        if ~isempty(idx)
            ocObj=varargin{idx};
        else
            error('No model is provided');
        end
        varargin(idx)=[];
        idx=find(cellfun('isclass',varargin,'double'));
        if length(idx)==3
            if numel(varargin{1})==1
                arcarg=varargin{1};
                varargin(1)=[];
            elseif numel(varargin{2}==1)
                arcarg=varargin{2};
                varargin(2)=[];
            elseif numel(varargin{3}==1)
                arcarg=varargin{3};
                varargin(3)=[];
            else
                error('')
            end
            if prod(size(varargin{1}))==numel(varargin{1})
                y=varargin{1};
                varargin(1)=[];
            elseif prod(size(varargin{2}))==numel(varargin{2})
                y=varargin{2};
                varargin(2)=[];
            else
                error('')
            end
            J=varargin{1};
            ocTrj.arcarg=arcarg;
            ocTrj.y=y(:);
            ocTrj.linearization=J;
        else
            error('')
        end
        ocTrj.modelname=modelname(ocObj);
        ocTrj.modelparameter=parametervalue(ocObj);
        dynPrim=dynprimitive(ocTrj);
end

function [ocTrj,period]=adddefault(dynPrimStruct)

if isfield(dynPrimStruct,'octrajectory') && isfield(dynPrimStruct,'period')
    ocTrj=dynPrimStruct.octrajectory;
    period=dynPrimStruct.period;
    return
elseif isfield(dynPrimStruct,'octrajectory')
    error('Field ''period'' is missing.')
end
mandatoryfields={'y','arcarg'};
optionalfields={'x','x0','linearization','solver','timehorizon','modelname','modelparameter','solverinfo','userinfo','violationinfo','arcposition','arcinterval'};
fn=fieldnames(dynPrimStruct);
for ii=1:length(mandatoryfields)
    if isempty(any(strcmp(mandatoryfields{ii},fn)))
        error('Field: is missing.\n')
    end
end
ocTrj.x=0;
ocTrj.arcinterval=[0 1];
ocTrj.arcposition=[1;1];
ocTrj.timehorizon=inf;
period=0;
for ii=1:length(fn)
    switch fn{ii}
        case optionalfields
            ocTrj.(fn{ii})=dynPrimStruct.(fn{ii});
        case mandatoryfields
            ocTrj.(fn{ii})=dynPrimStruct.(fn{ii});
    end
end
if size(ocTrj.y,2)~=1
    error('Input structure seems to be no equilibrium.')
end
ocTrj=octrajectory(ocTrj);