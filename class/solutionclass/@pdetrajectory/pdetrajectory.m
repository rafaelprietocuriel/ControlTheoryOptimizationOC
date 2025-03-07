function pdeTrj=pdetrajectory(varargin)
%
% PPDETRAJECTORY ppdetrajectory constructor
%
%
switch nargin
    case 0
        % generate empty pdeTrj
        pdeTrj.x=[];
        pdeTrj.y=[];
        pdeTrj.femdata=[];
        pdeTrj.arcarg=[];
        pdeTrj.arcposition=[];
        pdeTrj.arcinterval=[];

        pdeTrj.x0=[];
        pdeTrj.linearization=[];
        pdeTrj.solver='';
        pdeTrj.timehorizon=[];
        pdeTrj.modelname='';
        pdeTrj.modelparameter=[];
        pdeTrj.solverinfo=[];
        pdeTrj.userinfo=[];
        pdeTrj.violationinfo=[];
        pdeTrj=class(pdeTrj,'pdetrajectory');
        
    case 1
        % generate pdeTrj from struct or return pdetrajectory or empty
        % pdeTrj
        if ispdetrajectory(varargin{1})
            if ispdeprimitive(varargin{1}) || ispdeasymptotic(varargin{1})
                pdeStruct=struct(varargin{1});
                pdeTrj=pdeStruct.pdetrajectory;
            else
                pdeTrj=varargin{1};
            end
        elseif isstruct(varargin{1})
            [pdeTrj,mandatoryfields,optionalfields]=initfields(varargin{1});
            if isempty(pdeTrj)
                pdeTrj=pdetrajectory();
            end
            pdeTrj=orderfields(pdeTrj,[mandatoryfields optionalfields]);
            pdeTrj=class(pdeTrj,'pdetrajectory');
        elseif isempty(varargin{1})
            pdeTrj=pdetrajectory();
        else
            ocmaterror('Input argument is not a ppdetrajectory.')
        end
    case 2
        % add model information to pdeTrj, pdetrajectory(pdeTrj,pdeObj)
        if ispdetrajectory(varargin{1}) && isppdemodel(varargin{2})
            varargin{1}.modelname=modelname(varargin{2});
            varargin{1}.modelparameter=parametervalue(varargin{2});
            pdeTrj=pdetrajectory(varargin{1});
        elseif isstruct(varargin{1}) && isppdemodel(varargin{2})
            pdeTrj=pdetrajectory(varargin{1});
            if isempty(pdeTrj)
                return
            end
            pdeTrj=pdetrajectory(pdeTrj,varargin{2});
        end
end

function [pdeTrj,mandatoryfields,optionalfields]=initfields(pdeTrjInit)

% optional fields
optionalfields={'x0','linearization','solver','timehorizon','modelname','modelparameter','solverinfo','userinfo','violationinfo'};
mandatoryfields={'x','y','femdata','arcarg','arcposition','arcinterval'};

pdeTrj.x0=0;
pdeTrj.linearization=[];
pdeTrj.solver='';
pdeTrj.arcposition=[];
pdeTrj.timehorizon=[];
pdeTrj.modelname='';
pdeTrj.modelparameter=[];
pdeTrj.solverinfo=[];
pdeTrj.userinfo=[];
pdeTrj.violationinfo=[];

if ~all(ismember(mandatoryfields,fieldnames(pdeTrjInit)))
    idx=find(~ismember(mandatoryfields,fieldnames(pdeTrjInit)));
    ocmatmsg('Mandatory field ''%s'' is missing.\n',mandatoryfields{idx})
    pdeTrj=[];
    return
end
% mandatory fields
pdeTrj.x=pdeTrjInit.x;
pdeTrj.y=pdeTrjInit.y;
pdeTrj.femdata=pdeTrjInit.femdata;
pdeTrj.arcarg=pdeTrjInit.arcarg;
if ~isfield(pdeTrjInit,'arcposition')
    arcposition=find(diff(pdeTrjInit.x)==0);
    pdeTrjInit.arcposition=[1 arcposition+1;arcposition length(pdeTrjInit.x)];
end
pdeTrj.arcposition=pdeTrjInit.arcposition;
if isfield(pdeTrjInit,'solverinfo')
    pdeTrj.solverinfo=pdeTrjInit.solverinfo;
    pdeTrjInit=rmfield(pdeTrjInit,'solverinfo');
end
pdeTrj.arcinterval=pdeTrjInit.arcinterval;
if isfield(pdeTrjInit,'solver')
    pdeTrj.solver=pdeTrjInit.solver;
    pdeTrjInit=rmfield(pdeTrjInit,'solver');
end

pdeTrjInit=rmfield(pdeTrjInit,mandatoryfields);
% fieldnames that are neither of category mandatory or optional are put
% into field solverinfo
if ~isempty(pdeTrjInit)
    fn=fieldnames(pdeTrjInit);
    for ii=1:length(fn)
        pdeTrj.solverinfo.(fn{ii})=pdeTrjInit.(fn{ii});
    end
end