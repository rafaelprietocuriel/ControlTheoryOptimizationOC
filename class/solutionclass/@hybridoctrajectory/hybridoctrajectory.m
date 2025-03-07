function hocTrj=hybridoctrajectory(varargin)
%
% hybridoctrajectory hybridoctrajectory constructor
%
% hybridoctrajectory(ODESTRUCT) creates an hybridoctrajectory object from an ode solution
% structure ODESTRUCT.
%
% arcposition: starts at the second and the last but one entry of the
% fields x and y
% the first and last entries represent the left side and right side entries
% at the initial and end time
% jumparg: 0 ... no jump, 
%          1,2, ... jump with optimal jumptime
%          -1, -2, ...   jump with fixed jumptime 
%          -1.5, -2.5, ...   jump with fixed jumptime used for continuation
%                            from no jump to full jump, this is an
%                            intermediate  case and should be replaced by
%                            -1,-2, ... for the end trajectory.  
% 
switch nargin
    case 0
        hocTrj.x=[];
        hocTrj.y=[];
        hocTrj.arcarg=[];
        hocTrj.jumparg=[];
        hocTrj.arcposition=[];
        hocTrj.arcinterval=[];

        hocTrj.x0=[];
        hocTrj.linearization=[];
        hocTrj.solver='';
        hocTrj.timehorizon=[];
        hocTrj.modelname='';
        hocTrj.modelparameter=[];
        hocTrj.solverinfo=[];
        hocTrj.userinfo=[];
        hocTrj.violationinfo=[];
        hocTrj=class(hocTrj,'hybridoctrajectory');

    case 1
        if isstruct(varargin{1})
            % input argument is a structure consisting at least of the
            % mandatory fields
            try
                [hocTrj,mandatoryfields,optionalfields]=initfields(varargin{1});
                % remove mandatory fields
                name=fieldnames(varargin{1});
                for ii=1:numel(name) % if values for the optional fields are provided the default values are overwritten
                    if ismember(name{ii},optionalfields)
                        hocTrj.(name{ii})=varargin{1}.(name{ii});
                    else
                        hocTrj.solverinfo.(name{ii})=varargin{1}.(name{ii});
                    end
                end
                hocTrj=orderfields(hocTrj,[mandatoryfields optionalfields]);
                hocTrj=class(hocTrj,'hybridoctrajectory');
            catch
                rethrow(lasterror)
            end
        elseif isempty(varargin{1})
            hocTrj=hybridoctrajectory();
        elseif ishybridoctrajectory(varargin{1})
            hocTrj=hybridoctrajectory(struct(varargin{1}).hybridoctrajectory);
        elseif ishybridoctrajectory(varargin{1})
            if isempty(varargin{1})
                hocTrj=hybridoctrajectory();
            else
                hocTrj=varargin{1};
            end
        end
    case 2
        if isimpulsemodel(varargin{2})
            varargin{1}.modelname=modelname(varargin{2});
            varargin{1}.modelparameter=parametervalue(varargin{2});
            hocTrj=hybridoctrajectory(varargin{1});
        else
            ocmatmsg('Second argument is not an ocmodel and ignored.\n')
            hocTrj=hybridoctrajectory(varargin{1});
        end
end

function [hocTrj,mandatoryfields,optionalfields]=initfields(ocTrjInit)

% optional fields
hocTrj.x0=0;
hocTrj.linearization=[];
hocTrj.solver='';
hocTrj.timehorizon=[];
hocTrj.modelname='';
hocTrj.modelparameter=[];
hocTrj.solverinfo=[];
hocTrj.userinfo=[];
hocTrj.violationinfo=[];

optionalfields={'x0','linearization','solver','timehorizon','modelname','modelparameter','solverinfo','userinfo','violationinfo'};
% mandatory fields
hocTrj.x=ocTrjInit.x;
hocTrj.y=ocTrjInit.y;
hocTrj.arcarg=ocTrjInit.arcarg;
hocTrj.jumparg=ocTrjInit.jumparg;
if ~isfield(ocTrjInit,'arcposition')
    arcposition=find(diff(ocTrjInit.x)==0);
    ocTrjInit.arcposition=[1 arcposition+1;arcposition length(ocTrjInit.x)];
end
hocTrj.arcposition=ocTrjInit.arcposition;
if isfield(ocTrjInit,'solverinfo')
    hocTrj.solverinfo=ocTrjInit.solverinfo;
    ocTrjInit=rmfield(ocTrjInit,'solverinfo');
end
hocTrj.arcinterval=ocTrjInit.arcinterval;
if isfield(ocTrjInit,'solver')
    hocTrj.solver=ocTrjInit.solver;
    ocTrjInit=rmfield(ocTrjInit,'solver');
end

mandatoryfields={'x','y','arcarg','jumparg','arcposition','arcinterval'};
ocTrjInit=rmfield(ocTrjInit,mandatoryfields);
if ~isempty(ocTrjInit)
    fn=fieldnames(ocTrjInit);
    for ii=1:length(fn)
        hocTrj.solverinfo.(fn{ii})=ocTrjInit.(fn{ii});
    end
end