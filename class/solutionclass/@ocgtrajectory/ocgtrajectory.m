function ocgTrj=ocgtrajectory(varargin)
%
% OCGTRAJECTORY ocgtrajectory constructor
%
% OCGTRAJECTORY(ODESTRUCT) creates an octrajectory object from an ode solution
% structure ODESTRUCT.
%
idx=find(cellfun('isempty',varargin));
if ~isempty(idx)
    varargin(idx)=[];
    ocgTrj=ocgtrajectory(varargin{:});
    return
end
switch nargin
    case 0
        ocgTrj.odenum=[];
        ocgTrj=class(ocgTrj,'ocgtrajectory',octrajectory());
    case 1
        if isa(varargin{1},'ocgasymptotic')
            ocAsymS=struct(varargin{1});
            ocgTrj=ocAsymS.ocgtrajectory;
        elseif isa(varargin{1},'ocgtrajectory')
            ocgTrj=varargin{1};
        elseif isa(varargin{1},'ocasymptotic')
            ocTrj=octrajectory(varargin{1});
            arcpos=arcposition(ocTrj);
            n=arcnum(ocTrj);
            ocTrj=struct(ocTrj);
            y=ocTrj.y;
            odenum=size(y,1);
            ocTrj.y=[];
            for ii=1:n
                ocTrj.y{ii}=y(:,arcpos(1,ii):arcpos(2,ii));
            end
            ocTrj=octrajectory(ocTrj);
            ocgTrj.odenum=odenum(ones(1,n));
            ocgTrj=class(ocgTrj,'ocgtrajectory',ocTrj);
        elseif isa(varargin{1},'octrajectory')
            ocTrj=varargin{1};
            arcpos=arcposition(ocTrj);
            n=arcnum(ocTrj);
            ocTrj=struct(ocTrj);
            y=ocTrj.y;
            if isfield(ocTrj,'linearization')
                linearization=ocTrj.linearization;
            else
                linearization=[];
            end
            odenum=size(y,1);
            ocTrj.y=[];
            ocTrj.linearization=[];
            for ii=1:n
                ocTrj.y{ii}=y(:,arcpos(1,ii):arcpos(2,ii));
                if ~isempty(linearization)
                    if size(linearization,3)==1
                        ocTrj.linearization{ii}=linearization;
                    else
                    end
                end
            end
            ocTrj=octrajectory(ocTrj);
            ocgTrj.odenum=odenum(ones(1,n));
            ocgTrj=class(ocgTrj,'ocgtrajectory',ocTrj);
        elseif isstruct(varargin{1})
            if isfield(varargin{1},'odenum')
                odenum=varargin{1}.odenum;
                socTrj=rmfield(varargin{1},'odenum');
                ocTrj=octrajectory(socTrj);
                ocgTrj=ocgtrajectory(odenum,ocTrj);
            else
                ocgTrj=ocgtrajectory(octrajectory(varargin{1}));
            end
        else
            error('Input argument has to be an ''octrajectory'' a ''structure'' or an ''ocgtrajectory''.')
        end
    case 2
        
        idx=find(cellfun('isclass',varargin,'double'));
        if isempty(idx)
            error('No ''numeric''.')
        end
        odenum=varargin{idx};
        varargin(idx)=[];

        switch class(varargin{1})
            case 'octrajectory'
                ocTrj=varargin{1};
                odenum=odenum(:).';
                n=arcnum(ocTrj);
                arcpos=arcposition(ocTrj);
                if length(odenum)==n
                    if ~iscell(ocTrj.y)
                        y=ocTrj.y;
                        ocTrj.y=[];
                        for ii=1:n
                            ocTrj.y{ii}=y(1:odenum(ii),arcpos(1,ii):arcpos(2,ii));
                        end
                    end
                    ocTrj=octrajectory(ocTrj);
                    ocgTrj.odenum=odenum;
                    ocgTrj=class(ocgTrj,'ocgtrajectory',ocTrj);
                else
                    error('Wrong size of ''odenum''.')
                end
            case 'struct'
                ocTrj=octrajectory(varargin{1});
                ocgTrj=ocgtrajectory(ocTrj,odenum);
            otherwise
                error('Input arguments are not an ''octrajectory'' and ''numeric''.')

        end

    case 3
        idx=find(cellfun('isclass',varargin,'octrajectory'));
        if isempty(idx)
            idx=find(cellfun('isclass',varargin,'struct'));
            if isempty(idx)
                error('No ''octrajectory'' or ''struct''.')
            end
        end
        ocTrj=varargin{idx};
        varargin(idx)=[];
        idx=find(cellfun('isclass',varargin,'stdocmodel'));
        if isempty(idx)
            error('No ''stdocmodel''.')
        end
        ocObj=varargin{idx};
        varargin(idx)=[];
        odenum=varargin{1};
        ocTrj=octrajectory(ocTrj,ocObj);
        ocgTrj=ocgtrajectory(ocTrj,odenum);
end