function hocMultiPath=hybridocmultipath(varargin)
%
% OCMULTIPATH hybridocmultipath constructor
%
% OCMULTIPATH(ODESTRUCT) creates an hybridocmultipath object from an ode solution
% structure ODESTRUCT.
%
switch nargin
    case 0
        hocMultiPath.solutionclass{1}=hybridoctrajectory();
        hocMultiPath.solver='';
        hocMultiPath.solverinfo.tmesh=[];
        hocMultiPath.solverinfo.coeff=[];
        hocMultiPath.solverinfo.tangent=[];
        hocMultiPath.userinfo=[];
        hocMultiPath=class(hocMultiPath,'hybridocmultipath');
    case 1
        if isempty(varargin{1})
            hocMultiPath=hybridocmultipath();
        elseif ishybridoctrajectory(varargin{1})
            hocMultiPath=hybridocmultipath(varargin(1));
        elseif ishybridocmultipath(varargin{1})
            if isempty(varargin{1})
                hocMultiPath=hybridocmultipath();
            else
                hocMultiPath=varargin{1};
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
                    if ishybridoctrajectory(varargin{1}{ii})
                        if ~isempty(x0)
                            if abs(x0-initialtime(varargin{1}{ii}))>1e-4
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
                        ocmatmsg(['%d''th argument is not an hybridoctrajectory'])
                    end
                    solutionclass{multcounter}=varargin{1}{ii};
                end
            end
        end
        hocMultiPath.solutionclass=solutionclass;
        hocMultiPath.solver='';
        hocMultiPath.solverinfo.tmesh=[];
        hocMultiPath.solverinfo.coeff=[];
        hocMultiPath.solverinfo.tangent=[];
        if ~isempty(pathtype)
            hocMultiPath.solverinfo.pathtype=pathtype;
        end
        if ~isempty(inftimetransformation)
            hocMultiPath.solverinfo.inftimetransformation=inftimetransformation;
        end
        if ~isempty(multiarccalc)
            multiarccalc=unique(multiarccalc);
            if numel(multiarccalc)==1
                hocMultiPath.solverinfo.multiarccalc=multiarccalc;
            end
        end
        hocMultiPath.userinfo=[];
        hocMultiPath=class(hocMultiPath,'hybridocmultipath');
    otherwise
        hocMultiPathCell=[];
        multcounter=0;
        for ii=1:nargin
            if ~isempty(varargin{ii})
                multcounter=multcounter+1;
                if ishybridoctrajectory(varargin{ii})
                    hocMultiPathCell{multcounter}=varargin{ii};
                else
                    ocmatmsg(['%d''th argument is not an hybridoctrajectory'])
                end
            end
        end
        hocMultiPath=hybridocmultipath(hocMultiPathCell);
end