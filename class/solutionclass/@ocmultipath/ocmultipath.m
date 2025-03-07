function ocMultiPath=ocmultipath(varargin)
%
% OCMULTIPATH ocmultipath constructor
%
% OCMULTIPATH(ODESTRUCT) creates an ocmultipath object from an ode solution
% structure ODESTRUCT.
%
switch nargin
    case 0
        ocMultiPath.solutionclass{1}=octrajectory();
        ocMultiPath.solver='';
        ocMultiPath.solverinfo.tmesh=[];
        ocMultiPath.solverinfo.coeff=[];
        ocMultiPath.solverinfo.tangent=[];
        ocMultiPath.modelname='';
        ocMultiPath.modelparameter=[];
        ocMultiPath.userinfo=[];
        ocMultiPath=class(ocMultiPath,'ocmultipath');
    case 1
        if isempty(varargin{1})
            ocMultiPath=ocmultipath();
        elseif isocasymptotic(varargin{1}) || isoctrajectory(varargin{1})
            ocMultiPath=ocmultipath(varargin(1));
        elseif isocmultipath(varargin{1})
            if isempty(varargin{1})
                ocMultiPath=ocmultipath();
            else
                ocMultiPath=varargin{1};
                return
            end
        elseif iscell(varargin{1})
            multcounter=0;
            % mandatoryfields
            x0=[];
            modelpar=[];
            pathtype=[];
            inftimetransformation=[];
            multiarccalc=[];
            % remove empty cells
            varargin{1}(cellfun('isempty',varargin{1}))=[];
            numvarargin=numel(varargin{1});
            solutionclass=cell(1,numvarargin);
            modelnam=cell(1,numvarargin);
            modelpar=cell(1,numvarargin);
            for ii=1:numvarargin
                if ~isempty(varargin{1}{ii})
                    multcounter=multcounter+1;
                    if isocasymptotic(varargin{1}{ii}) || isoctrajectory(varargin{1}{ii})
                        if ~isempty(x0)
                            if abs(x0-initialtime(varargin{1}{ii}))>1e-4
                                ocmaterror('Initial times are not consistent.')
                            end
                        else
                            x0=initialtime(varargin{1}{ii});
                        end
                        modelnam{ii}=modelname(varargin{1}{ii});
                        modelpar{ii}=modelparameter(varargin{1}{ii});
                        solverinf=solverinfo(varargin{1}{ii});
                        if isfield(solverinf,'pathtype')
                            pathtype{ii}=solverinf.pathtype;
                        end
                        if isfield(solverinf,'inftimetransformation')
                            if isempty(solverinf.inftimetransformation)
                                solverinf.inftimetransformation=0;
                            end
                            inftimetransformation(ii)=solverinf.inftimetransformation;
                        end
                        if isfield(solverinf,'multiarccalc')
                            multiarccalc(ii)=solverinf.multiarccalc;
                        end
                    else
                        ocmatmsg(['%d''th argument is not an ocasymptotic or octrajectory'])
                    end
                    solutionclass{multcounter}=varargin{1}{ii};
                end
            end
        end
        ocMultiPath.solutionclass=solutionclass;
        ocMultiPath.solver='';
        ocMultiPath.solverinfo.tmesh=[];
        ocMultiPath.solverinfo.coeff=[];
        ocMultiPath.solverinfo.tangent=[];
        if ~isempty(pathtype)
            ocMultiPath.solverinfo.pathtype=pathtype;
        end
        ocMultiPath.modelname=modelnam;
        ocMultiPath.modelparameter=modelpar;
        if ~isempty(inftimetransformation)
            ocMultiPath.solverinfo.inftimetransformation=inftimetransformation;
        end
        if ~isempty(multiarccalc)
            multiarccalc=unique(multiarccalc);
            if numel(multiarccalc)==1
                ocMultiPath.solverinfo.multiarccalc=multiarccalc;
            end
        end
        ocMultiPath.userinfo=[];
        ocMultiPath=class(ocMultiPath,'ocmultipath');
    otherwise
        ocMultiPathCell=[];
        multcounter=0;
        for ii=1:nargin
            if ~isempty(varargin{ii})
                multcounter=multcounter+1;
                if isocasymptotic(varargin{ii}) || isoctrajectory(varargin{ii})
                    ocMultiPathCell{multcounter}=varargin{ii};
                else
                    ocmatmsg(['%d''th argument is not an ocasymptotic or octrajectory'])
                end
            end
        end
        ocMultiPath=ocmultipath(ocMultiPathCell);
end