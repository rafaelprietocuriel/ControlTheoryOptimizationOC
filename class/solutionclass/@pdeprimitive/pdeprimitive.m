function pdePrim=pdeprimitive(varargin)
%
%

switch nargin
    case 0
        pdeTrj=pdetrajectory([]);
        pdePrim.period=[];
        pdePrim=class(pdePrim,'pdeprimitive',pdeTrj);
        
    case 1
        if isempty(varargin{1})
            pdePrim=pdeprimitive();
        elseif ispdeprimitive(varargin{1})
            pdePrim=varargin{1};
        elseif isstruct(varargin{1})
            if isfield(varargin{1},'pdetrajectory')
                pdeTrj=varargin{1}.pdetrajectory;
            else
                if isfield(varargin{1},'period')
                    varargin{1}=rmfield(varargin{1},'period');
                    pdeTrj=varargin{1};
                else
                    pdeTrj=varargin{1};
                end
            end
            if isstruct(pdeTrj)
                pdeTrj=adddefault(pdeTrj);
                pdeTrj=pdetrajectory(pdeTrj);
            end
            pdePrim=pdeprimitive(pdeTrj);
        elseif ispdeasymptotic(varargin{1})
            pdePrim=varargin{1}.pdeprimitive;
        elseif ispdetrajectory(varargin{1})
            T=arcinterval(varargin{1});
            if length(varargin{1}.x)==1
                T=0;
            else
                if isinf(T(end))
                    T=0;
                else
                    T=T(end);
                end
            end
            pdePrim.period=T;
            pdeTrj=pdetrajectory(varargin{1});
            pdePrim=class(pdePrim,'pdeprimitive',pdeTrj);
        end
    case 2
        if isppdemodel(varargin{2}) && ispdetrajectory(varargin{1})
            varargin{1}.modelname=modelname(varargin{2});
            varargin{1}.modelparameter=parametervalue(varargin{2});
            pdePrim=pdeprimitive(varargin{1});
        elseif isstruct(varargin{1}) && isppdemodel(varargin{2})
            varargin{1}.modelname=modelname(varargin{2});
            varargin{1}.modelparameter=parametervalue(varargin{2});
            pdePrim=pdeprimitive(varargin{1});
        end
        
    case 3
        if isnumeric(varargin{1}) && isnumeric(varargin{2}) && numel(varargin{2})==1 && isstruct(varargin{3})
            % equilibrium given by its values and arc argument
            pdePrimStruct.y=varargin{1};
            pdePrimStruct.arcarg=varargin{2};
            pdePrimStruct.femdata=varargin{3};
            pdePrim=pdeprimitive(pdePrimStruct);
        end
        
end


function pdeTrjStruct=adddefault(pdeTrjStruct)

if isfield(pdeTrjStruct,'pdetrajectory') % for equilibrium case default values can be used
    pdeTrjStruct=pdeTrjStruct.pdetrajectory;
    pdeTrjStruct=rmfield(pdeTrjStruct,'pdetrajectory');
end
if size(pdeTrjStruct.y,2)==1
    pdeTrjStruct.x=0;
    pdeTrjStruct.arcposition=[1;1];
    pdeTrjStruct.arcinterval=[0 1];
    pdeTrjStruct.timehorizon=inf;
end
