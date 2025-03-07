function gradCuv=gradcurve(varargin)
%
% GRADCURVE gradcurve constructor
%

switch nargin
    case 0
        gradCuv.continuationtype='';
        gradCuv.continuationparameter=[];
        gradCuv.objectivevalue=[];
        gradCuv.continuationdata=[];

        gradCuv.modelname='';
        gradCuv.modelparameter=[];
        gradCuv.curveinfo=[];
        gradCuv.userinfo=[];
        gradCuv=class(gradCuv,'gradcurve');

    case 1
        if isstruct(varargin{1})
            try
                [gradCuv,optionalfields]=initoptionalfields();
                gradCuv.objectivevalue=varargin{1}.objectivevalue;
                gradCuv.continuationtype=varargin{1}.continuationtype;
                gradCuv.continuationdata=varargin{1}.continuationdata;
                gradCuv.continuationparameter=varargin{1}.continuationparameter;
                varargin{1}=rmfield(varargin{1},{'objectivevalue','continuationparameter','continuationdata','continuationtype'});
                name=fieldnames(varargin{1});
                for ii=1:numel(name)
                    if ismember(name{ii},optionalfields)
                        gradCuv.(name{ii})=varargin{1}.(name{ii});
                    end
                end
                gradCuv=orderfields(gradCuv,[{'continuationtype','continuationparameter','objectivevalue','continuationdata'} optionalfields]);
                gradCuv=class(gradCuv,'gradcurve');
            catch
                ocmatmsg('\nStructure of ''%s'' is not consistent with structure of an ''gradcurve''.\n\n',inputname(1))
                lasterr
            end
        elseif isempty(varargin{1})
            gradCuv=gradcurve();
        elseif isoccurve(varargin{1})
            if isempty(varargin{1})
                gradCuv=gradcurve();
            else
                gradCuv=varargin{1};
            end
        end
    case 2
        if isocmodel(varargin{2})
            varargin{1}.modelname=modelname(varargin{2});
            varargin{1}.modelparameter=parametervalue(varargin{2});
            gradCuv=gradcurve(varargin{1});
        else
            ocmatmsg('Second argument is not an ocmodel and ignored.\n')
            gradCuv=gradcurve(varargin{1});
        end
end

function [gradCuv,optionalfields]=initoptionalfields()

gradCuv.modelname='';
gradCuv.modelparameter=[];
gradCuv.curveinfo=[];
gradCuv.userinfo=[];

optionalfields=fieldnames(gradCuv).';