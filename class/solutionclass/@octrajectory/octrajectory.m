function ocTrj=octrajectory(varargin)
%
% OCTRAJECTORY octrajectory constructor
%
% OCTRAJECTORY(ODESTRUCT) creates an octrajectory object from an ode solution
% structure ODESTRUCT.
%
idx=find(cellfun('isempty',varargin));
if ~isempty(idx)
    varargin(idx)=[];
    ocTrj=octrajectory(varargin{:});
    return
end

switch nargin
    case 0
        ocTrj=class(initfields(),'octrajectory');

    case 1
        if isstruct(varargin{1})
            [ocTrj,mandatoryfields,optionalfields]=initfields();
            fn=fieldnames(varargin{1});
            for ii=1:length(mandatoryfields)
                if ~(any(strcmp(mandatoryfields{ii},fn)))
                    if strcmp(mandatoryfields{ii},'arcposition')
                        arcposition=find(diff(varargin{1}.x)==0);
                        varargin{1}.arcposition=[1 arcposition+1; ...
                            arcposition length(varargin{1}.x)];
                    else
                        ocmaterror('Field: %s is missing in octrjactory structure.\n',mandatoryfields{ii})
                    end
                end
            end
            fn=fieldnames(varargin{1});
            for ii=1:length(fn)
                switch fn{ii}
                    case optionalfields
                        ocTrj.(fn{ii})=varargin{1}.(fn{ii});
                    case mandatoryfields
                        ocTrj.(fn{ii})=varargin{1}.(fn{ii});
                end
            end
            if isempty(ocTrj.x0)
                ocTrj.x0=0;
            end
            ocTrj=orderfields(ocTrj,[mandatoryfields optionalfields]);
            ocTrj=class(ocTrj,'octrajectory');
        elseif isdynprimitive(varargin{1}) ||  isocasymptotic(varargin{1})
            ocTrj=octrajectory(struct(varargin{1}).octrajectory);
        elseif isoctrajectory(varargin{1})
            ocTrj=varargin{1};
        end
    case 2
        idx=find(cellfun('isclass',varargin,'stdocmodel'));
        if ~isempty(idx)
            ocObj=varargin{idx};
        else
            ocObj=[];
        end
        varargin(idx)=[];
        idx=find(cellfun('isclass',varargin,'double'));
        if ~isempty(idx)
            J=varargin{idx};
        else
            J=[];
        end
        varargin(idx)=[];
        if isempty(varargin)
            error('No ''octrajectory provided''')
        end
        
        if ~isempty(ocObj)
            varargin{1}.modelname=modelname(ocObj);
            varargin{1}.modelparameter=parametervalue(ocObj);
            ocTrj=octrajectory(varargin{1});
        elseif ~isempty(J)
            ocTrj=varargin{1};
            ocTrj.linearization=J;
        else
            ocmatmsg('Type of second argument is unknown and ignored.\n')
            ocTrj=octrajectory(varargin{1});
        end
end

function [ocTrj,mandatoryfields,optionalfields]=initfields()

ocTrj.x=[];
ocTrj.y=[];
ocTrj.arcarg=[];
ocTrj.arcposition=[];
ocTrj.arcinterval=[];

ocTrj.x0=[];
ocTrj.linearization=[];
ocTrj.solver='';
ocTrj.timehorizon=[];
ocTrj.modelname='';
ocTrj.modelparameter=[];
ocTrj.solverinfo=[];
ocTrj.userinfo=[];
ocTrj.violationinfo=[];

optionalfields={'x0','linearization','solver','timehorizon','modelname','modelparameter','solverinfo','userinfo','violationinfo'};
mandatoryfields={'x','y','arcarg','arcposition','arcinterval'};
