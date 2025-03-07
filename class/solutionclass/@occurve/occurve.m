function ocCuv=occurve(varargin)
%
% OCCURVE occurve constructor
%

switch nargin
    case 0
        ocCuv.y=[];
        ocCuv.arcarg=[];
        ocCuv.arcposition=[];

        ocCuv.modelname='';
        ocCuv.modelparameter=[];
        ocCuv.curveinfo=[];
        ocCuv.userinfo=[];
        ocCuv=class(ocCuv,'occurve');

    case 1
        if isstruct(varargin{1})
            try
                [ocCuv,optionalfields]=initoptionalfields();
                ocCuv.y=varargin{1}.y;
                ocCuv.arcarg=varargin{1}.arcarg;
                ocCuv.arcposition=varargin{1}.arcposition;
                varargin{1}=rmfield(varargin{1},{'y','arcarg','arcposition'});
                name=fieldnames(varargin{1});
                for ii=1:numel(name)
                    if ismember(name{ii},optionalfields)
                        ocCuv.(name{ii})=varargin{1}.(name{ii});
                    end
                end
                ocCuv=orderfields(ocCuv,[{'y','arcarg','arcposition'} optionalfields]);
                ocCuv=class(ocCuv,'occurve');
            catch
                ocmatmsg('\nStructure of ''%s'' is not consistent with structure of an ''occurve''.\n\n',inputname(1))
                lasterr
            end
        elseif isempty(varargin{1})
            ocCuv=occurve();
        elseif isnumeric(varargin{1})
            ocCuv.y=varargin{1};
            ocCuv.arcarg=0;
            ocCuv.arcposition=[1 size(varargin{1},2)]';
            ocCuv=occurve(ocCuv);
        elseif isoccurve(varargin{1})
            if isempty(varargin{1})
                ocCuv=occurve();
            else
                ocCuv=varargin{1};
            end
        elseif isocasymptotic(varargin{1})
            ocCuv=occurve(octrajectory(varargin{1}));
        elseif isoctrajectory(varargin{1})
            ocCuv=occurve(struct(varargin{1}));
        elseif isdynprimitive(varargin{1})
            ocCuv=occurve(octrajectory(varargin{1}));
        elseif iscell(varargin{1}) && numel(varargin{1})==1
            if isoccurve(varargin{1}{1})
                ocCuv=varargin{1}{1};
            end
        end
    case 2
        if isocmodel(varargin{2})
            varargin{1}.modelname=modelname(varargin{2});
            varargin{1}.modelparameter=parametervalue(varargin{2});
            ocCuv=occurve(varargin{1});
        else
            ocmatmsg('Second argument is not an ocmodel and ignored.\n')
            ocCuv=occurve(varargin{1});
        end
end

function [ocCuv,optionalfields]=initoptionalfields()

ocCuv.modelname='';
ocCuv.modelparameter=[];
ocCuv.curveinfo=[];
ocCuv.userinfo=[];

optionalfields=fieldnames(ocCuv).';