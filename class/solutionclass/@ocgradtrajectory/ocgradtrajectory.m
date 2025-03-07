function ocgTrj=ocgradtrajectory(varargin)
%
% OCGRADTRAJECTORY ocgradtrajectory constructor
%
% OCGRADTRAJECTORY(ODESTRUCT) creates an ocgradtrajectory object from an ode solution
% structure ODESTRUCT.
%
switch nargin
    case 0
        ocgTrj.variable=struct(); % y,z,q,cst_y,cst_z,cst_q,u,v,w,o
        ocgTrj.argument=struct(); % t ... time, a ... age
        ocgTrj.initargument=struct(); % t0, a0
        ocgTrj.objectivevalue=[];
        ocgTrj.tangent=struct();
        ocgTrj.type=''; % possible types are: loc ... local (non-distributed) variables, age ... age distributed, space ... spatial distributed
        ocgTrj.solver=struct();
        ocgTrj.timehorizon=[];
        ocgTrj.modelname='';
        ocgTrj.modelparameter=[];
        ocgTrj.solverinfo=struct();
        ocgTrj.userinfo=struct();
        ocgTrj.violationinfo=struct();
        ocgTrj=class(ocgTrj,'ocgradtrajectory');
        
    case 1
        if isstruct(varargin{1})
            if isfield(varargin{1},'argument')
                % input argument is a structure consisting at least of the
                % mandatory fields
                try
                    [ocgTrj,mandatoryfields,optionalfields]=initfields(varargin{1});
                    % remove mandatory fields
                    name=fieldnames(varargin{1});
                    for ii=1:numel(name) % if values for the optional fields are provided the default values are overwritten
                        if ismember(name{ii},optionalfields)
                            if ~isempty(varargin{1}.(name{ii}))
                                ocgTrj.(name{ii})=varargin{1}.(name{ii});
                            end
                        else
                            ocgTrj.solverinfo.(name{ii})=varargin{1}.(name{ii});
                        end
                    end
                    ocgTrj.type=determineoctype(ocgTrj);
                    ocgTrj=orderfields(ocgTrj,[mandatoryfields optionalfields]);
                    ocgTrj=class(ocgTrj,'ocgradtrajectory');
                catch
                    rethrow(lasterror)
                end
            elseif isfield(varargin{1},'extremal')
                ocgTrj=ocgradtrajectory(varargin{1}.extremal,varargin{2:end});
            else
                [ocgTrj,mandatoryfields,optionalfields,argumentfields,initargumentfields,variablefields]=initfields(varargin{1});
                name=fieldnames(varargin{1});
                for ii=1:numel(name) % if values for the optional fields are provided the default values are overwritten
                    if ismember(name{ii},argumentfields)
                        ocgTrj.argument.(name{ii})=varargin{1}.(name{ii});
                    elseif ismember(name{ii},variablefields)
                        ocgTrj.variable.(name{ii})=varargin{1}.(name{ii});
                    elseif ismember(name{ii},initargumentfields)
                        ocgTrj.initargument.(name{ii})=varargin{1}.(name{ii});
                    elseif ismember(name{ii},optionalfields)
                        ocgTrj.(name{ii})=varargin{1}.(name{ii});
                    else
                        ocgTrj.solverinfo.(name{ii})=varargin{1}.(name{ii});
                    end
                end
                fn=fieldnames(ocgTrj);
                ocgTrj=orderfields(ocgTrj,[mandatoryfields optionalfields]);
                ocgTrj=ocgradtrajectory(ocgTrj);
            end
        elseif isempty(varargin{1})
            ocgTrj=ocgradtrajectory();
        elseif isocgradtrajectory(varargin{1})
                ocgTrj=varargin{1};
        elseif ischar(varargin{1})
            ocgTrj=ocgradtrajectory();
            switch varargin{1}
                case 'concentrated'
                    variable.y=[];
                    variable.cst_y=[];
                    variable.v=[];
                    variable.o=[];
                    argument.t=[];
                    initargument.t0=[];
                    ocgTrj.variable=variable;
                    ocgTrj.argument=argument;
                    ocgTrj.initargument=initargument;
                    ocgTrj.type=varargin{1};
                case 'age'
                case 'space'
            end
        end
end

function [ocgTrj,mandatoryfields,optionalfields,argumentfields,initargumentfields,variablefields]=initfields(ocTrjInit)
mandatoryfields={'variable'};

% optional fields
ocgTrj.objectivevalue=[];
ocgTrj.tangent=struct();
ocgTrj.type='';
ocgTrj.solver=struct();
ocgTrj.timehorizon=[];
ocgTrj.argument=[];
ocgTrj.initargument=[];
ocgTrj.modelname='';
ocgTrj.modelparameter=[];
ocgTrj.solverinfo=struct();
ocgTrj.userinfo=struct();
ocgTrj.violationinfo=[];

optionalfields={'argument','initargument','objectivevalue','tangent','type','solver','timehorizon','modelname','modelparameter','solverinfo','userinfo','violationinfo'};
if isfield(ocTrjInit,'initargument')
    % mandatory fields
    ocgTrj.initargument=ocTrjInit.initargument;
    ocgTrj.argument=ocTrjInit.argument;
    ocgTrj.variable=ocTrjInit.variable;
    if isfield(ocTrjInit,'solverinfo')
        ocgTrj.solverinfo=ocTrjInit.solverinfo;
        ocTrjInit=rmfield(ocTrjInit,'solverinfo');
    end
    ocTrjInit=rmfield(ocTrjInit,mandatoryfields);
    if ~isempty(ocTrjInit)
        fn=fieldnames(ocTrjInit);
        for ii=1:length(fn)
            ocgTrj.solverinfo.(fn{ii})=ocTrjInit.(fn{ii});
        end
    end
    if ~isempty(ocgTrj.argument)
        ocgTrj.timehorizon=ocgTrj.argument.t(end);
    end
end
variablefields={'u','v','w','z','y','q','cst_z','cst_y','cst_q','o'};
argumentfields={'t','a'};
initargumentfields={'t0','a0'};

function octype=determineoctype(ocgTrj)

if isempty(ocgTrj.argument)
    octype='static';
    return
end
name=fieldnames(ocgTrj.variable);
flag=[1 1];
for ii=1:numel(name) % if values for the optional fields are provided the default values are overwritten
    if ismember(name{ii},{'u','w','z','q'})
        flag(2)=3;
    elseif ismember(name{ii},{'v','y'})
        flag(1)=2;
    end
end
if prod(flag)==2
    octype='concentrated';
elseif prod(flag)==3 || prod(flag)==6
    octype='age';
end
    


