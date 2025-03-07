function ocgMTrj=ocgradmultipath(varargin)
%
% OCGRADTRAJECTORY ocgradtrajectory constructor
%
% OCGRADTRAJECTORY(ODESTRUCT) creates an ocgradtrajectory object from an ode solution
% structure ODESTRUCT.
%
switch nargin
    case 0
        ocgMTrj.argument=struct(); % t ... time, a ... age
        ocgMTrj.variable=struct(); % y,z,q,cst_y,cst_z,cst_q,u,v,w,o
        ocgMTrj.initargument=struct(); % t0, a0
        ocgMTrj.multiplicity=struct(); % multiplicity.number ... number of solutions, multiplicity.timestructure=[1 a2 ...;b1 b2 ...] ... first row starting indices for the different solutions
        ocgMTrj.objectivevalue=[];
        ocgMTrj.tangent=struct();
        ocgMTrj.type=''; % possible types are: loc ... local (non-distributed) variables, age ... age distributed, space ... spatial distributed
        ocgMTrj.solver=struct();
        ocgMTrj.timehorizon=[];
        ocgMTrj.modelname='';
        ocgMTrj.modelparameter=[];
        ocgMTrj.solverinfo=struct();
        ocgMTrj.userinfo=struct();
        ocgMTrj.violationinfo=struct();
        ocgMTrj=class(ocgMTrj,'ocgradmultipath');

    case 1
        if isempty(varargin{1})
            ocgMTrj=ocgradmultipath();
        elseif isocgradmultipath(varargin{1})
            ocgMTrj=varargin{1};
        elseif isocgradtrajectory(varargin{1})
            if isempty(varargin{1})
                ocgMTrj=ocgradmultipath();
                return
            end
            ocgMTrj=struct(varargin{1});
            ocgMTrj.multiplicity.number=1;
            switch octype(varargin{1})
                case 'concentrated'
                    ocgMTrj.multiplicity.timestructure=[1;length(ocgMTrj.argument.t)];
            end
            ocgMTrj=orderfields(ocgMTrj,orderedfieldnames());
            ocgMTrj=class(ocgMTrj,'ocgradmultipath');
        elseif iscell(varargin{1})
            ocgMTrj=ocgradmultipath(varargin{1}{:});
        elseif isstruct(varargin{1})
            if length(varargin{1})>1
                ocgTrjC=cell(1,length(varargin{1}));
                for ii=1:length(varargin{1});
                    ocgTrjC{ii}=ocgradtrajectory(varargin{1}(ii));
                end
                ocgMTrj=ocgradmultipath(ocgTrjC{:});
            else
                ocgMTrj=ocgradtrajectory(varargin{1});
            end
        end
    otherwise
        ctr=0;
        for ii=1:nargin
            if isempty(varargin{ii})
                ocgMTrj=ocgradmultipath();
                ocmatmsg('Return empty ''multipath'' since one argument is empty.')
                return
            end
            ctr=ctr+1;
            if isocgradtrajectory(varargin{ii})
                if ii==1
                    ocgMTrj=struct(varargin{ii});
                    ocgMTrj.multiplicity.number=ctr;
                    switch octype(varargin{ii})
                        case 'concentrated'
                            ocgMTrj.multiplicity.timestructure=[1;length(ocgMTrj.argument.t)];
                        case 'static'
                            ocgMTrj.multiplicity.timestructure=[];
                    end
                    act_modelname=modelname(varargin{ii});
                    act_modelparameter=modelparameter(varargin{ii});
                else
                    ocgTrj=struct(varargin{ii});
                    if ~strcmp(act_modelname,modelname(varargin{ii}))
                        ocgMTrj=ocgradmultipath();
                        ocmatmsg('Return ''ocgradtrajectory'' are solutions of different models.')
                        return
                    end
                    if length(act_modelparameter)~=length(modelparameter(varargin{ii}))% || any(act_modelparameter-modelparameter(varargin{ii}))
                        ocgMTrj=ocgradmultipath();
                        ocmatmsg('Return ''ocgradtrajectory'' are computed for different parameter values.')
                        return
                    end
                    ocgMTrj.multiplicity.number=ctr;
                    if strcmp(ocgMTrj.type,'static')
                        ocgMTrj.argument=[];
                        ocgMTrj.initargument=[];
                        ocgMTrj.timehorizon=[];
                    else
                        ocgMTrj.multiplicity.timestructure(:,ctr)=ocgMTrj.multiplicity.timestructure(2,ctr-1)+[1;length(ocgTrj.argument.t)];
                        ocgMTrj.argument.t=[ocgMTrj.argument.t ocgTrj.argument.t];
                        ocgMTrj.initargument.t0=[ocgMTrj.initargument.t0 ocgTrj.initargument.t0];
                        ocgMTrj.timehorizon(ctr)=ocgTrj.timehorizon;
                    end
                    fn=fieldnames(ocgTrj.variable);
                    for jj=1:length(fn)
                        ocgMTrj.variable.(fn{jj})=[ocgMTrj.variable.(fn{jj}) ocgTrj.variable.(fn{jj})];
                    end
                    ocgMTrj.objectivevalue(ctr)=ocgTrj.objectivevalue;
                end
            else
                ocgMTrj=ocgradmultipath();
                ocmatmsg('Return empty ''multipath'' since one argument is not a ''ocgradtrajectory''.')
                return
            end
        end
        ocgMTrj=orderfields(ocgMTrj,orderedfieldnames());
        ocgMTrj=class(ocgMTrj,'ocgradmultipath');
end

function fn=orderedfieldnames()
fn=fieldnames(struct(ocgradmultipath()));
function [ocgMTrj,mandatoryfields,optionalfields,argumentfields,initargumentfields,variablefields]=initfields(ocTrjInit)
mandatoryfields={'argument','variable','initargument'};

% optional fields
ocgMTrj.objectivevalue=[];
ocgMTrj.tangent=struct();
ocgMTrj.type='';
ocgMTrj.solver=struct();
ocgMTrj.timehorizon=[];
ocgMTrj.modelname='';
ocgMTrj.modelparameter=[];
ocgMTrj.solverinfo=struct();
ocgMTrj.userinfo=struct();
ocgMTrj.violationinfo=[];

optionalfields={'objectivevalue','tangent','type','solver','timehorizon','modelname','modelparameter','solverinfo','userinfo','violationinfo'};
if isfield(ocTrjInit,'initargument')
    % mandatory fields
    ocgMTrj.initargument=ocTrjInit.initargument;
    ocgMTrj.argument=ocTrjInit.argument;
    ocgMTrj.variable=ocTrjInit.variable;
    if isfield(ocTrjInit,'solverinfo')
        ocgMTrj.solverinfo=ocTrjInit.solverinfo;
        ocTrjInit=rmfield(ocTrjInit,'solverinfo');
    end
    ocTrjInit=rmfield(ocTrjInit,mandatoryfields);
    if ~isempty(ocTrjInit)
        fn=fieldnames(ocTrjInit);
        for ii=1:length(fn)
            ocgMTrj.solverinfo.(fn{ii})=ocTrjInit.(fn{ii});
        end
    end
    if ~isempty(ocgMTrj)
        ocgMTrj.timehorizon=ocgMTrj.argument.t(end);
    end
end
variablefields={'u','v','w','z','y','q','cst_z','cst_y','cst_q','o'};
argumentfields={'t','a'};
initargumentfields={'t0','a0'};

function octype=determineoctype(ocgMTrj)

name=fieldnames(ocgMTrj.variable);
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
    


