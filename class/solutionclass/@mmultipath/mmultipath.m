function ocMultiPath=mmultipath(varargin)
%
% MMULTIPATH mmultipath constructor
%
% MMULTIPATH(ODESTRUCT) creates an mmultipath object from an ode solution
% structure ODESTRUCT.
%
% USERINFO: has to contain the information if mmultipath is an indifference
% solution. The fields are MULTIPLICITY, ARCSTRUCTURE
switch nargin
    case 0
        ocMultiPath.solutionclass{1}=octrajectory();
        ocMultiPath.parts=[];
        ocMultiPath.optimaltransition=[];
        ocMultiPath.solver='';
        ocMultiPath.solverinfo=struct([]);
        ocMultiPath.modelname=cell('');
        ocMultiPath.modelparameter=[];
        ocMultiPath.userinfo=[];
        ocMultiPath=class(ocMultiPath,'mmultipath');
    case 1
        if isempty(varargin{1})
            ocMultiPath=mmultipath();
        elseif isocasymptotic(varargin{1}) || isoctrajectory(varargin{1})
            ocMultiPath=mmultipath(varargin(1));
        elseif ismmultipath(varargin{1})
            if isempty(varargin{1})
                ocMultiPath=mmultipath();
            else
                ocMultiPath=varargin{1};
                return
            end
        elseif iscell(varargin{1})
            varargin{1}(cellfun('isempty',varargin{1}))=[];
            solutionparts=numel(varargin{1});
            solutionclass=cell(1,solutionparts);
            modelnam=cell(1,solutionparts);
            modelpar=cell(1,solutionparts);
            optimaltransition=zeros(1,solutionparts);
            solverinf=cell(1,solutionparts);
            ctr=0;
            for ii=1:solutionparts
                if isocasymptotic(varargin{1}{ii}) || isoctrajectory(varargin{1}{ii})
                    ctr=ctr+1;
                    solutionclass{ctr}=varargin{1}{ii};
                    modelnam{ctr}=modelname(varargin{1}{ii});
                    modelpar{ctr}=modelparameter(varargin{1}{ii});
                    solverinf{ctr}=solverinfo(varargin{1}{ii});
                    optimaltransition(ctr)=0;
                elseif ismmultipath(varargin{1}{ii})
                    ctr0=ctr+1;
                    ctr=ctr+length(varargin{1}{ii}.solutionclass);
                    solutionclass(ctr0:ctr)=varargin{1}{ii}.solutionclass;
                    modelnam(ctr0:ctr)=varargin{1}{ii}.modelname;
                    modelpar(ctr0:ctr)=varargin{1}{ii}.modelparameter;
                    solverinf(ctr0:ctr)=varargin{1}{ii}.solverinfo;
                    optimaltransition(ctr0:ctr)=varargin{1}{ii}.optimaltransition;
                else
                    
                    ocmatmsg(['%d''th argument is not an ocasymptotic or octrajectory'])
                end
            end
            ocMultiPath.solutionclass=solutionclass;
            ocMultiPath.parts=ctr;
            ocMultiPath.optimaltransition=optimaltransition;
            ocMultiPath.solver='';
            ocMultiPath.solverinfo=solverinf;
            ocMultiPath.modelname=modelnam;
            ocMultiPath.modelparameter=modelpar;
            ocMultiPath.userinfo=[];
            ocMultiPath=class(ocMultiPath,'mmultipath');
        end
    otherwise
        if iscell(varargin{1}) || ismmultipath(varargin{1})
            ocMultiPath=mmultipath(varargin{1});
            for ii=2:2:nargin
                switch varargin{ii}
                    case 'optimaltransition'
                        ocMultiPath.optimaltransition=varargin{ii+1};
                    case 'userinfo'
                        if isempty(ocMultiPath.userinfo)
                            ocMultiPath.userinfo=varargin{ii+1};
                        else
                            fn=fieldnames(varargin{ii+1});
                            for jj=1:length(fn)
                                ocMultiPath.userinfo.(fn{jj})=varargin{ii+1}.(fn{jj});
                            end
                        end
                    otherwise
                    ocmatmsg(['input has to be either ''userinfo'' or ''optimaltransition''.'])
                end
            end
        else
        end
end