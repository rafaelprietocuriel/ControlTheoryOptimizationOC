function mapPrim=mapprimitive(varargin)
%
%

switch nargin
    case 0
        mapPrim.period=[];
        docTrj=doctrajectory([]);
        mapPrim=class(mapPrim,'mapprimitive',docTrj);

    case 1
        if isstruct(varargin{1})
            try
                varargin{1}=adddefault(varargin{1});
                docTrj=doctrajectory(varargin{1}.doctrajectory);
                %varargin{1}=rmfield(varargin{1},'doctrajectory');
                if ~isfield(varargin{1},'period') && size(docTrj.y,2)==1
                    varargin{1}.period=0;
                end
                mapPrim.period=varargin{1}.period;
                mapPrim=class(mapPrim,'mapprimitive',docTrj);
            catch
                rethrow(lasterror)
            end
        elseif ismapprimitive(varargin{1})
            if isempty(varargin{1})
                mapPrim=mapprimitive();
            else
                mapPrim=varargin{1};
            end
        elseif isempty(varargin{1})
            mapPrim=mapprimitive();
        end
    case 2
        if (isocmodel(varargin{2}) || isdocmodel(varargin{2})) || isodemodel(varargin{2}) && isstruct(varargin{1})
            varargin{1}=adddefault(varargin{1});
            varargin{1}.doctrajectory=doctrajectory(varargin{1}.doctrajectory,varargin{2});
            mapPrim=mapprimitive(varargin{1});
        elseif isnumeric(varargin{1}) && isnumeric(varargin{2}) && numel(varargin{2})==1
            % equilibrium given by its values and arc argument
            mapPrimStruct.y=varargin{1};
            mapPrimStruct.arcarg=varargin{2};
            mapPrim=mapprimitive(mapPrimStruct);
        elseif  isnumeric(varargin{1}) && isoctrajectory(varargin{2})
            mapPrimStruct.period=varargin{1};
            mapPrimStruct.doctrajectory=doctrajectory(varargin{2});
            mapPrim=mapprimitive(mapPrimStruct);
        else
            ocmatmsg('Second argument is not an ocmodel and ignored.\n')
            mapPrim=mapprimitive(varargin{1});
        end
    case 3

        if isnumeric(varargin{1}) && isnumeric(varargin{2}) ...
                && numel(varargin{2})==1 && isocmodel(varargin{3})
            mapPrimStruct.y=varargin{1};
            mapPrimStruct.arcarg=varargin{2};
            mapPrim=mapprimitive(mapPrimStruct,varargin{3});
        elseif isnumeric(varargin{1}) && isnumeric(varargin{2}) ...
                && numel(varargin{2})==1 && isnumeric(varargin{3})
            mapPrimStruct.y=varargin{1};
            mapPrimStruct.arcarg=varargin{2};
            mapPrimStruct.linearization=varargin{3};
            mapPrim=mapprimitive(mapPrimStruct);
        else
        end
end

function mapPrimStruct=adddefault(mapPrimStruct)

if ~isfield(mapPrimStruct,'doctrajectory') % for equilibrium case default values can be used
    if size(mapPrimStruct.y,2)==1
        mapPrimStruct.doctrajectory.x=0;
        mapPrimStruct.doctrajectory.y=mapPrimStruct.y;
        mapPrimStruct.doctrajectory.x0=0;
        mapPrimStruct.doctrajectory.y0=mapPrimStruct.y;
        mapPrimStruct.doctrajectory.arcarg=mapPrimStruct.arcarg;
        mapPrimStruct.doctrajectory.arcposition=[1 1]';
        mapPrimStruct.doctrajectory.timehorizon=inf;
        mapPrimStruct.period=1;
        if isfield(mapPrimStruct,'linearization')
            mapPrimStruct.doctrajectory.linearization=mapPrimStruct.linearization;
        end
    else
        if ~isfield(mapPrimStruct,'doctrajectory')
            try
                mapPrimStruct.doctrajectory.x=mapPrimStruct.x;
                mapPrimStruct.doctrajectory.y=mapPrimStruct.y;
                mapPrimStruct.doctrajectory.x0=mapPrimStruct.x0;
                mapPrimStruct.doctrajectory.y0=mapPrimStruct.y0;
                mapPrimStruct.doctrajectory.arcarg=mapPrimStruct.arcarg;
                mapPrimStruct.doctrajectory.arcposition=[1 size(mapPrimStruct.y,2)]';
                mapPrimStruct.doctrajectory.timehorizon=inf;
                if isfield(mapPrimStruct,'linearization')
                    mapPrimStruct.doctrajectory.linearization=mapPrimStruct.linearization;
                end
            end
        end
    end
end