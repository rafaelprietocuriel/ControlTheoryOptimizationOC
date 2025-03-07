function ppdeTrj=ppdetrajectory(varargin)
%
% PPDETRAJECTORY ppdetrajectory constructor
%
% PPDETRAJECTORY(ODESTRUCT) creates an ppdetrajectory object from an ode solution
% structure ODESTRUCT.
%
switch nargin
    case 0
        ppdeTrj.femdata=[];
        %ppdeTrj.x=[];
        ppdeTrj=class(ppdeTrj,'ppdetrajectory',octrajectory());
        
    case 1
        if isppdetrajectory(varargin{1})
            ppdeTrj=varargin{1};
        elseif isstruct(varargin{1})
            if isfield(varargin{1},'femdata')
                ocTrj=varargin{1};
                ocTrj=rmfield(ocTrj,'femdata');
                ppdeTrj=ppdetrajectory(ocTrj,varargin{1}.femdata);
            else
            end
        elseif isempty(varargin{1})
            ppdeTrj=ppdetrajectory();
        else
            ocmaterror('Input argument is not a ppdetrajectory.')
        end
    case 2
        if isppdemodel(varargin{2})
            varargin{1}.modelname=modelname(varargin{2});
            varargin{1}.modelparameter=parametervalue(varargin{2});
            ppdeTrj=ppdetrajectory(varargin{1});
        elseif isoctrajectory(varargin{1}) && isstruct(varargin{2})
            if isempty(varargin{1})
                ppdeTrj=ppdetrajectory();
            else
                ppdeTrj.femdata=varargin{2};
                ppdeTrj=class(ppdeTrj,'ppdetrajectory',varargin{1});
            end
        elseif isoctrajectory(varargin{2}) && isstruct(varargin{1})
            ppdeTrj.femdata=varargin{1};
            ppdeTrj=class(ppdeTrj,'ppdetrajectory',varargin{2});
        elseif isstruct(varargin{1}) && isstruct(varargin{2})
            try
                ocTrj=octrajectory(varargin{1});
                ppdeTrj.femdata=varargin{2};
            catch
                try
                    ocTrj=octrajectory(varargin{2});
                    ppdeTrj.femdata=varargin{1};
                catch
                    ocmaterror('Input argument is not an octrajectory and a femdata structure.');
                end
            end
            ppdeTrj=class(ppdeTrj,'ppdetrajectory',ocTrj);
        elseif isppdetrajectory(varargin{1}) && isnumeric(varargin{2})
            ppdeTrj=varargin{1};
            if size(dependentvar(ppdeTrj),1)==size(varargin{2},1) && size(varargin{2},1)==size(varargin{2},2)
                ppdeTrj.octrajectory=octrajectory(ppdeTrj.octrajectory,varargin{2});
            end
        end
end

function [ocTrj,mandatoryfields,optionalfields]=initfields(ocTrjInit)

% optional fields
ocTrj.x0=0;
ocTrj.linearization=[];
ocTrj.solver='';
ocTrj.timehorizon=[];
ocTrj.modelname='';
ocTrj.modelparameter=[];
ocTrj.solverinfo=[];
ocTrj.userinfo=[];
ocTrj.violationinfo=[];

optionalfields={'x0','linearization','solver','timehorizon','modelname','modelparameter','solverinfo','userinfo','violationinfo'};
% mandatory fields
ocTrj.x=ocTrjInit.x;
ocTrj.y=ocTrjInit.y;
ocTrj.arcarg=ocTrjInit.arcarg;
if ~isfield(ocTrjInit,'arcposition')
    arcposition=find(diff(ocTrjInit.x)==0);
    ocTrjInit.arcposition=[1 arcposition+1;arcposition length(ocTrjInit.x)];
end
ocTrj.arcposition=ocTrjInit.arcposition;
if isfield(ocTrjInit,'solverinfo')
    ocTrj.solverinfo=ocTrjInit.solverinfo;
    ocTrjInit=rmfield(ocTrjInit,'solverinfo');
end
ocTrj.arcinterval=ocTrjInit.arcinterval;
if isfield(ocTrjInit,'solver')
    ocTrj.solver=ocTrjInit.solver;
    ocTrjInit=rmfield(ocTrjInit,'solver');
end

mandatoryfields={'x','y','arcarg','arcposition','arcinterval'};
ocTrjInit=rmfield(ocTrjInit,mandatoryfields);
if ~isempty(ocTrjInit)
    fn=fieldnames(ocTrjInit);
    for ii=1:length(fn)
        ocTrj.solverinfo.(fn{ii})=ocTrjInit.(fn{ii});
    end
end