function ppdePrim=ppdeprimitive(varargin)
%
%

switch nargin
    case 0
        ppdePrim.period=[];
        ocTrj=ppdetrajectory([]);
        ppdePrim=class(ppdePrim,'ppdeprimitive',ocTrj);

    case 1
        if isstruct(varargin{1})
            try
                if ~isfield(varargin{1},'femdata')
                    ppdePrim=ppdeprimitive();
                    return
                end
                varargin{1}=adddefault(varargin{1});
                ocTrj=ppdetrajectory(varargin{1}.ppdetrajectory);
                %varargin{1}=rmfield(varargin{1},'ppdetrajectory');
                if ~isfield(varargin{1},'period')
                    if size(dependentvar(ocTrj),2)==1
                        varargin{1}.period=0;
                    else
                        varargin{1}.period=ocTrj.arcinterval(end);
                    end
                end
                ppdePrim.period=varargin{1}.period;
                ppdePrim=class(ppdePrim,'ppdeprimitive',ocTrj);
            catch
                rethrow(lasterror)
            end
        elseif isppdeprimitive(varargin{1})
            if isempty(varargin{1})
                ppdePrim=ppdeprimitive();
            else
                ppdePrim=varargin{1};
            end
        elseif isocasymptotic(varargin{1})
            ppdePrimStruct.ppdetrajectory=struct(ppdetrajectory(varargin{1}));
            ppdePrimStruct.ppdetrajectory.solverinfo=[];
            ppdePrimStruct.period=ppdePrimStruct.ppdetrajectory.arcinterval(end);
            ppdePrim=ppdeprimitive(ppdePrimStruct);
        elseif isppdetrajectory(varargin{1})
            ppdePrimStruct=struct(varargin{1});
            ppdePrimStruct.period=ppdePrimStruct.arcinterval(end);
            ppdePrim=ppdeprimitive(ppdePrimStruct);
        elseif isempty(varargin{1})
            ppdePrim=ppdeprimitive();
        end
    case 2
        if isppdemodel(varargin{2}) && isstruct(varargin{1})
            varargin{1}=adddefault(varargin{1});
            varargin{1}.ppdetrajectory=ppdetrajectory(varargin{1}.ppdetrajectory,varargin{2});
            ppdePrim=ppdeprimitive(varargin{1});
        elseif isnumeric(varargin{1}) && isnumeric(varargin{2}) && numel(varargin{2})==1
            % equilibrium given by its values and arc argument
            ppdePrimStruct.y=varargin{1};
            ppdePrimStruct.arcarg=varargin{2};
            ppdePrim=ppdeprimitive(ppdePrimStruct);
        elseif isppdeprimitive(varargin{1}) && isnumeric(varargin{2})
            ppdePrim=varargin{1};
            if period(varargin{1})==0
                if size(dependentvar(ppdePrim),1)==size(varargin{2},1) && size(varargin{2},1)==size(varargin{2},2)
                    ppdePrim.ppdetrajectory=ppdetrajectory(ppdePrim.ppdetrajectory,varargin{2});
                end
            else
            end
        end
    case 3

        if isnumeric(varargin{1}) && isnumeric(varargin{2}) ...
                && numel(varargin{2})==1 && isocmodel(varargin{3})
            ppdePrimStruct.y=varargin{1};
            ppdePrimStruct.arcarg=varargin{2};
            ppdePrim=ppdeprimitive(ppdePrimStruct,varargin{3});
        elseif isnumeric(varargin{1}) && isnumeric(varargin{2}) ...
                && numel(varargin{2})==1 && isnumeric(varargin{3})
            ppdePrimStruct.y=varargin{1};
            ppdePrimStruct.arcarg=varargin{2};
            ppdePrimStruct.linearization=varargin{3};
            ppdePrim=ppdeprimitive(ppdePrimStruct);
        else
        end
end

function ppdePrimStruct=adddefault(ppdePrimStruct)

if ~isfield(ppdePrimStruct,'ppdetrajectory') % for equilibrium case default values can be used
    if size(ppdePrimStruct.y,2)==1
        ppdePrimStruct.ppdetrajectory.x=0;
        ppdePrimStruct.ppdetrajectory.y=ppdePrimStruct.y;
        ppdePrimStruct.ppdetrajectory.arcarg=ppdePrimStruct.arcarg;
        ppdePrimStruct.ppdetrajectory.femdata=ppdePrimStruct.femdata;
        ppdePrimStruct.ppdetrajectory.arcinterval=[0 1];
        ppdePrimStruct.ppdetrajectory.timehorizon=inf;
        ppdePrimStruct.period=0;
        if isfield(ppdePrimStruct,'linearization')
            ppdePrimStruct.ppdetrajectory.linearization=ppdePrimStruct.linearization;
        end
        if isfield(ppdePrimStruct,'userinfo')
            ppdePrimStruct.ppdetrajectory.userinfo=ppdePrimStruct.userinfo;
        end
    end
end