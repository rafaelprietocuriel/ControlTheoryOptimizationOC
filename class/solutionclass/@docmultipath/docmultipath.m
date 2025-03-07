function ocMultiPath=docmultipath(varargin)
%
% OCMULTIPATH ocmultipath constructor
%
% OCMULTIPATH(ODESTRUCT) creates an ocmultipath object from an ode solution
% structure ODESTRUCT.
%
switch nargin
    case 0
        ocMultiPath.solutionclass{1}=doctrajectory();
        ocMultiPath.solver='';
        ocMultiPath.solverinfo.tmesh=[];
        ocMultiPath.solverinfo.coeff=[];
        ocMultiPath.solverinfo.tangent=[];
        ocMultiPath.userinfo=[];
        ocMultiPath=class(ocMultiPath,'docmultipath');
    case 1
        if isempty(varargin{1})
            ocMultiPath=docmultipath();
        elseif isdocasymptotic(varargin{1}) || isdoctrajectory(varargin{1})
            ocMultiPath=docmultipath(varargin(1));
        elseif isdocmultipath(varargin{1})
            if isempty(varargin{1})
                ocMultiPath=docmultipath();
            else
                ocMultiPath=varargin{1};
                return
            end
        elseif iscell(varargin{1})
            multcounter=0;
            % mandatoryfields
            x0=[];
            modelnam='';
            modelpar=[];
            pathtype=[];
            inftimetransformation=[];
            multiarccalc=[];
            % remove empty cells
            varargin{1}(cellfun('isempty',varargin{1}))=[];
            numvarargin=numel(varargin{1});
            solutionclass=cell(1,numvarargin);
            for ii=1:numvarargin
                if ~isempty(varargin{1}{ii})
                    multcounter=multcounter+1;
                    if isdocasymptotic(varargin{1}{ii}) || isdoctrajectory(varargin{1}{ii})
                        if ~isempty(x0)
                            if x0~=initialtime(varargin{1}{ii})
                                ocmaterror('Initial times are not consistent.')
                            end
                        else
                            x0=initialtime(varargin{1}{ii});
                        end
                        if ~isempty(modelnam)
                            if ~strcmp(modelnam,modelname(varargin{1}{ii}))
                                ocmaterror('Models are not consistent.')
                            end
                        else
                            modelnam=modelname(varargin{1}{ii});
                        end
                        if ~isempty(modelpar)
                            try
                                if ~all(modelpar==modelparameter(varargin{1}{ii}))
                                    ocmaterror('Parametervalues are not consistent.')
                                end
                            catch
                                ocmaterror('Parametervalues are not consistent.')
                            end
                        else
                            modelpar=modelparameter(varargin{1}{ii});
                        end
                        solverinf=solverinfo(varargin{1}{ii});
                        if isfield(solverinf,'pathtype')
                            pathtype{ii}=solverinf.pathtype;
                        end
                        if isfield(solverinf,'multiarccalc')
                            multiarccalc(ii)=solverinf.multiarccalc;
                        end
                    else
                        ocmatmsg(['%d''th argument is not an ocasymptotic or doctrajectory'])
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
        ocMultiPath=class(ocMultiPath,'docmultipath');
    otherwise
        ocMultiPathCell=[];
        multcounter=0;
        for ii=1:nargin
            if ~isempty(varargin{ii})
                multcounter=multcounter+1;
                if isdocasymptotic(varargin{ii}) || isdoctrajectory(varargin{ii})
                    ocMultiPathCell{multcounter}=varargin{ii};
                else
                    ocmatmsg(['%d''th argument is not an ocasymptotic or doctrajectory'])
                end
            end
        end
        ocMultiPath=docmultipath(ocMultiPathCell);
end