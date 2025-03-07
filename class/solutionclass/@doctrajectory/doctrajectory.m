function docTrj=doctrajectory(varargin)
%
% OCTRAJECTORY doctrajectory constructor
%
% OCTRAJECTORY(ODESTRUCT) creates an doctrajectory object from an ode solution
% structure ODESTRUCT.
%
switch nargin
    case 0
        docTrj.x=[];
        docTrj.y=[];
        docTrj.arcarg=[];
        docTrj.arcposition=[];

        docTrj.x0=[];
        docTrj.y0=[];
        docTrj.linearization=[];
        docTrj.solver='';
        docTrj.timehorizon=[];
        docTrj.modelname='';
        docTrj.modelparameter=[];
        docTrj.solverinfo=[];
        docTrj.userinfo=[];
        docTrj.violationinfo=[];
        docTrj=class(docTrj,'doctrajectory');

    case 1
        if isstruct(varargin{1})
            % input argument is a structure consisting at least of the
            % mandatory fields
            try
                [docTrj,mandatoryfields,optionalfields]=initfields(varargin{1});
                % remove mandatory fields
                name=fieldnames(varargin{1});
                for ii=1:numel(name) % if values for the optional fields are provided the default values are overwritten
                    if ismember(name{ii},optionalfields)
                        docTrj.(name{ii})=varargin{1}.(name{ii});
                    end
                end
                docTrj=orderfields(docTrj,[mandatoryfields optionalfields]);
                docTrj=class(docTrj,'doctrajectory');
            catch
                rethrow(lasterror)
            end
        elseif isempty(varargin{1})
            docTrj=doctrajectory();
        elseif ismapprimitive(varargin{1}) ||  isdocasymptotic(varargin{1})
            docTrj=doctrajectory(struct(varargin{1}).doctrajectory);
        elseif isdoctrajectory(varargin{1})
            if isempty(varargin{1})
                docTrj=doctrajectory();
            else
                docTrj=varargin{1};
            end
        end
    case 2
        if isdocmodel(varargin{2})
            varargin{1}.modelname=modelname(varargin{2});
            varargin{1}.modelparameter=parametervalue(varargin{2});
            docTrj=doctrajectory(varargin{1});
        else
            ocmatmsg('Second argument is not an ocmodel and ignored.\n')
            docTrj=doctrajectory(varargin{1});
        end
end

function [docTrj,mandatoryfields,optionalfields]=initfields(ocTrjInit)

% optional fields
docTrj.x0=0;
docTrj.y0=[];
docTrj.linearization=[];
docTrj.solver='';
docTrj.timehorizon=[];
docTrj.modelname='';
docTrj.modelparameter=[];
docTrj.solverinfo=[];
docTrj.userinfo=[];
docTrj.violationinfo=[];

optionalfields={'x0','y0','linearization','solver','timehorizon','modelname','modelparameter','solverinfo','userinfo','violationinfo'};
% mandatory fields
docTrj.x=ocTrjInit.x;
docTrj.y=ocTrjInit.y;
docTrj.arcarg=ocTrjInit.arcarg;
if ~isfield(ocTrjInit,'arcposition')
    arcposition=find(diff(ocTrjInit.x)==0);
    ocTrjInit.arcposition=[1 arcposition+1;arcposition length(ocTrjInit.x)];
end
docTrj.arcposition=ocTrjInit.arcposition;
if isfield(ocTrjInit,'solverinfo')
    docTrj.solverinfo=ocTrjInit.solverinfo;
    ocTrjInit=rmfield(ocTrjInit,'solverinfo');
end
mandatoryfields={'x','y','arcarg','arcposition'};
ocTrjInit=rmfield(ocTrjInit,mandatoryfields);
if ~isempty(ocTrjInit)
    fn=fieldnames(ocTrjInit);
    for ii=1:length(fn)
        docTrj.solverinfo.(fn{ii})=ocTrjInit.(fn{ii});
    end
end