function statEx=staticextremal(varargin)
%
% OCGRADTRAJECTORY ocgradtrajectory constructor
%
% OCGRADTRAJECTORY(ODESTRUCT) creates an ocgradtrajectory object from an ode solution
% structure ODESTRUCT.
%
switch nargin
    case 0
        statEx.x=[]; % y,z,q,cst_y,cst_z,cst_q,u,v,w,o
        statEx.lm=[]; % t ... time, a ... age
        statEx.activeindex=[];
        statEx.objectivevalue=[];
        statEx.tangent=struct();
        statEx.type=''; % possible types are: loc ... local (non-distributed) variables, age ... age distributed, space ... spatial distributed
        statEx.solver=struct();
        statEx.modelname='';
        statEx.modelparameter=[];
        statEx.solverinfo=struct();
        statEx.userinfo=struct();
        statEx.violationinfo=struct();
        statEx=class(statEx,'staticextremal');
        
    case 1
        if isstruct(varargin{1})
            if isfield(varargin{1},'x')
                % input argument is a structure consisting at least of the
                % mandatory fields
                try
                    [statEx,mandatoryfields,optionalfields]=initfields(varargin{1});
                    % remove mandatory fields
                    name=fieldnames(varargin{1});
                    for ii=1:numel(name) % if values for the optional fields are provided the default values are overwritten
                        if ismember(name{ii},optionalfields)
                            if ~isempty(varargin{1}.(name{ii}))
                                statEx.(name{ii})=varargin{1}.(name{ii});
                            end
                        else
                            statEx.solverinfo.(name{ii})=varargin{1}.(name{ii});
                        end
                    end
                    statEx=orderfields(statEx,[mandatoryfields optionalfields]);
                    statEx=class(statEx,'staticextremal');
                catch
                    rethrow(lasterror)
                end
            end
        elseif isempty(varargin{1})
            statEx=staticextremal();
        elseif isstaticextremal(varargin{1})
                statEx=varargin{1};
        end
end

function [statEx,mandatoryfields,optionalfields]=initfields(statExInit)
mandatoryfields={'x','lm'};
statEx.x=statExInit.x;
statEx.lm=statExInit.lm;
statExInit=rmfield(statExInit,mandatoryfields);

% optional fields
statEx.activeindex=[];
statEx.objectivevalue=[];
statEx.tangent=struct();
statEx.type=''; % possible types are: loc ... local (non-distributed) variables, age ... age distributed, space ... spatial distributed
statEx.solver=struct();
statEx.modelname='';
statEx.modelparameter=[];
statEx.solverinfo=struct();
statEx.userinfo=struct();
statEx.violationinfo=struct();

optionalfields={'activeindex','objectivevalue','tangent','type','solver','modelname','modelparameter','solverinfo','userinfo','violationinfo'};
for ii=1:length(optionalfields)
    if isfield(statExInit,optionalfields{ii})
        statEx.(optionalfields{ii})=statExInit.(optionalfields{ii});
        statExInit=rmfield(statExInit,optionalfields{ii});
    end
end
fn=fieldnames(statExInit);
for ii=1:length(fn)
        statEx.solverinfo.(fn{ii})=statExInit.(fn{ii});
end
