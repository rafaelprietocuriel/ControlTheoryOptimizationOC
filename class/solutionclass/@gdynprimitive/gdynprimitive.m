function gdynPrim=gdynprimitive(varargin)
%
%

idx=find(cellfun('isempty',varargin));
if ~isempty(idx)
    varargin(idx)=[];
    gdynPrim=gdynprimitive(varargin{:});
    return
end

switch nargin
    case 0
        gdynPrim.period=[];
        ocgTrj=ocgtrajectory([]);
        gdynPrim=class(gdynPrim,'gdynprimitive',ocgTrj);

    case 1
        switch class(varargin{1})
            case 'dynprimitive'
                gdynPrim.period=varargin{1}.period;
                ocgTrj=ocgtrajectory(varargin{1}.octrajectory);
                gdynPrim=class(gdynPrim,'gdynprimitive',ocgTrj);
            case 'gdynprimitive'
                gdynPrim=varargin{1};
            case 'struct'
                if isfield(varargin{1},'period')
                    period=varargin{1}.period;
                    varargin{1}=rmfield(varargin{1},'period');
                else
                    period=0;
                end
                if isfield(varargin{1},'odenum')
                    odenum=varargin{1}.odenum;
                    varargin{1}=rmfield(varargin{1},'odenum');
                else
                    odenum=[];
                end
                gdynPrim.period=period;
                ocgTrj=ocgtrajectory(varargin{1},odenum);
                gdynPrim=class(gdynPrim,'gdynprimitive',ocgTrj);
        end
    case 2
    case 3
end
