function varargout=store(dgObj,varargin)
%
% STORE results of continuation process
%
% STORE(dgObj0,PROBLEMTYPE) store results of continuation process
% given by PROBLEMTYPE:
%   'extremal2ep'   ... saddle-path of an equilibrium, continuing along the
%                       initial point
%   'extremalp2ep'  ... saddle-path of an equilibrium, continuing along a
%                       parameter value
%   'extremalt2ep'  ... saddle-path of an equilibrium, continuing the
%                       truncation time
%   'indifferencesolution'  ... continuation of an indifference threshold
%                       (Skiba point)
%   'limitextremal' ... continuation of a limitpoint solution
%
% STORE(dgObj0,OCELEMENT) stores the object OCELEMENT in the Result field of
% the stdocmodel object dgObj0. 
% Possible classes for OCELEMENT are
%        dynprimitive: equilibrium or limit cycle
%        ocasymptotic: asymptotic solutions
%        cell array of 'dynprimitives' or 'ocasymptotics'
%
% The default names of the Result field are
%       Result.Equilibrium: for a dynprimitive being an equilibrium
%       Result.LimitCycle: for a dynprimitive being a limit cycle
%       Result.ExtremalSolution: for an ocasymptotic
%       Result.Continuation: for the result of a continuation process
% If no output argument is specified it is written into the object dgObj0.
%
% dgObj0N=STORE(...) the result is returned to the new stdocmodel instance
% dgObj0N 


if isempty(dgObj) || nargin==1
    return
end

fieldname=[];
fieldvalue=[];
if ~ischar(varargin{1}) && nargin==2
    if isdynprimitive(varargin{1}) ||  isocasymptotic(varargin{1}) || isoctrajectory(varargin{1}) || isoccurve(varargin{1}) ||  isocmultipath(varargin{1})
        fieldname=class(varargin{1});
    elseif iscell(varargin{1})
        for ii=1:numel(varargin{1})
            store(dgObj,varargin{1}{ii},varargin{2:end});
        end

        if ~nargout
            assignin('caller',inputname(1),dgObj);
        else
            varargout{1}=dgObj;
        end
        return
    else
        ocmaterror('Input argument ''%s'' needs the specification of a field name.\n',inputname(2))
    end
    fieldvalue=varargin{1};
end

if ischar(varargin{1}) && nargin==2
    fieldname=varargin{1};
    fieldvalue=[];
end

if ischar(varargin{1}) && nargin==3
    fieldname=varargin{1};
    fieldvalue=varargin{2};
end

if isempty(fieldname) && isempty(fieldvalue)
    ocmaterror('Wrong input arguments')
end
[ocStruct,ocResultFieldName]=generateelement(dgObj,fieldname,fieldvalue);
if isempty(ocResultFieldName) || isempty(ocStruct)
    return
end
if isfield(dgObj.Result,ocResultFieldName)
    dgObj.Result.(ocResultFieldName){end+1}=ocStruct;
else
    dgObj.Result(1).(ocResultFieldName){1}=ocStruct;
end

if ~nargout
    assignin('caller',inputname(1),dgObj);
else
    varargout{1}=dgObj;
end

%%
function [ocResultStruct,ocResultFieldName]=generateelement(dgObj,fieldname,fieldvalue)
ocResultStruct=[];
ocResultFieldName='';
switch fieldname
    case 'dynprimitive'
        if isequilibrium(fieldvalue)
            ocResultStruct=fieldvalue;
            ocResultFieldName='Equilibrium';
        else
            ocResultStruct=fieldvalue;
            ocResultFieldName='PeriodicSolution';
        end
    case 'ocasymptotic'
        ocResultStruct=fieldvalue;
        ocResultFieldName='ExtremalSolution';
        
    case {'extremaldg2ep','extremal2ep'}
        global OCMATCONT OCMATAE OCBVP

        [globalvarfile,resultfile,OCMATCONTORIGINAL,OCMATAEORIGINAL,OCBVPORIGINAL]=loadglobal(dgObj,fieldname);

        % load data of information files into workspace
        load(globalvarfile)
        load(resultfile)

        mainfunch=feval(str2func(fieldname));

        % determine if non admissible solution occurs
        [soltarget,sollast,soln,sol]=findspecificsolution(sout);
        counter=0;
        for ii=1:numel(bvpout)
            try
                counter=counter+1;
                DataAdaptation(bvpout(ii).tmesh);
                ocExStruct(counter)=mainfunch{22}(bvpout(ii).tmesh,bvpout(ii).coeff,bvpout(ii).tangent);
                slMfStruct.y(:,counter)=ocExStruct(counter).y(:,1);
                if isfield(OCMATAE,'objectivevaluecoord') & ~isempty(OCMATAE.objectivevaluecoord)
                    slMfStruct.y(OCMATAE.objectivevaluecoord,counter)=ocExStruct(counter).y(OCMATAE.objectivevaluecoord,end);
                end
                slMfStruct.arcarg(counter)=OCMATCONT.HE.arcarg(1);
                slMfStruct.userinfo.continuationparameter(counter)=bvpout(ii).coeff(end);
            end
        end
        slMfStruct.arcposition=[];
        if ~isempty(OCMATAE.freeparameter)

        else
            ocEP=generateequilibrium(OCMATAE.EP.saddlepoint,OCMATAE.EP.arcarg,parametervalue(dgObj),modelname(dgObj), ...
                'dynprimitive',0,'differentialgame');
            if ~isempty(soln)
                ocEPN=ocEP;
            else
                ocEPN=[];
            end
        end
        if ~isempty(soltarget)
            sol=soltarget;
        else
            sol=ocExStruct(end);
        end
        if ~isempty(soln)
            ocResultStruct.NonadmissibleSolution=ocasymptotic(octrajectory(soln),ocEPN);
        else
            ocResultStruct.NonadmissibleSolution=ocasymptotic([]);
        end
        ocResultStruct.ExtremalSolution=ocasymptotic(octrajectory(sol),ocEP);

        ocResultStruct.SliceManifold=occurve(slMfStruct);
        ocResultStruct.ContinuationSolution=ocExStruct;
        ocResultStruct.ContinuationClassification=fieldname;
        ocResultStruct.ContinuationInformation=sout;
        ocResultStruct.LimitSet=ocEP;
        ocResultFieldName='Continuation';

        resetglobal(fieldname,OCMATCONTORIGINAL,OCMATAEORIGINAL,OCBVPORIGINAL);

    case {'indifferencesolution'}
        global OCMATCONT OCMATINDIF OCBVP

        [globalvarfile,resultfile,OCMATCONTORIGINAL,OCMATINDIFORIGINAL,OCBVPORIGINAL]=loadglobal(dgObj,fieldname);

        % load data of information files into workspace
        load(globalvarfile)
        load(resultfile)

        mainfunch=feval(str2func(fieldname));


        % determine if non admissible solution occurs
        [soltarget,sollast,soln,sol]=findspecificsolution(sout);
        counter=0;
        for ii=1:numel(bvpout)
            try
                counter=counter+1;
                DataAdaptation(bvpout(ii).tmesh);
                ocExStruct(counter)=mainfunch{22}(bvpout(ii).tmesh,bvpout(ii).coeff,bvpout(ii).tangent);
                indifferencepointStruct.y(1:OCMATCONT.maxnumode,counter)=ocExStruct(ii).y(1:OCMATCONT.maxnumode,1);
                indifferencepointStruct.arcarg(counter)=OCMATCONT.HE.arcarg(1);
                indifferencepointStruct.userinfo.continuationparameter(counter)=bvpout(ii).coeff(end);
            end
        end
        indifferencepointStruct.arcposition=[];

        %         ocEP=dynprimitive(OCMATINDIF.saddlepoint,OCMATCONT.HE.arcarg(end),OCMATINDIF.linearization);
        if ~isempty(soln)
            if any(strcmp(fieldname,{'indifferencesolution','indifferencedistribution','indifferencesolution4ae2ftae'}))
                switch soln.solver
                    case 'gbvp4c'
                        ocResultStruct.NonadmissibleSolution=sol2ocmultipath(soln,OCMATINDIF,OCMATINDIF.limitset,OCMATINDIF.linearization);
                        ocResultStruct.ExtremalSolution=sol2ocmultipath(sol,OCMATINDIF,OCMATINDIF.limitset,OCMATINDIF.linearization);
                    otherwise
                        ocResultStruct.NonadmissibleSolution=sol2ocmultipath(soln,OCMATINDIF,OCMATINDIF.saddlepoint,OCMATINDIF.linearization);
                        ocResultStruct.ExtremalSolution=sol2ocmultipath(sol,OCMATINDIF,OCMATINDIF.saddlepoint,OCMATINDIF.linearization);
                end
            elseif strcmp(fieldname,'indifferencesolution4per')
                ocResultStruct.NonadmissibleSolution=sol2ocmultipath(soln,OCMATINDIF,OCMATINDIF.limitset,OCMATINDIF.linearization);
                ocResultStruct.ExtremalSolution=sol2ocmultipath(sol,OCMATINDIF,OCMATINDIF.limitset,OCMATINDIF.linearization);
            elseif strcmp(fieldname,'indifferencesolution4emf')
                ocResultStruct.NonadmissibleSolution=sol2ocmultipath(soln,OCMATINDIF,soln.solverinfo.saddlepoint,OCMATINDIF.linearization);
                ocResultStruct.ExtremalSolution=sol2ocmultipath(sol,OCMATINDIF,sol.solverinfo.saddlepoint,OCMATINDIF.linearization);
            elseif strcmp(fieldname,'indifferencesolution4ft')
                ocResultStruct.NonadmissibleSolution=sol2ocmultipath(soln);
                ocResultStruct.ExtremalSolution=sol2ocmultipath(sol);
            elseif strcmp(fieldname,'indifferencesolutionep') 
                dgObj0tmp=changeparametervalue(dgObj0,OCMATINDIF.freeparameterindex,sol.solverinfo.parameters(OCMATINDIF.freeparametercoordinate));
                hatx.y=sol.solverinfo.parameters(sol.solverinfo.equilibriumcoord);
                hatx.x=0;
                hatx.arcarg=OCMATCONT.HE.arcarg(end);
                hatx.linearization=jacobian(dgObj0tmp,hatx);
                ocEP=dynprimitive(hatx);
                ocResultStruct.ExtremalSolution=ocasymptotic(octrajectory(sol,dgObj0tmp),ocEP);
                dgObj0tmp=changeparametervalue(dgObj0,OCMATINDIF.freeparameterindex,soln.solverinfo.parameters(OCMATINDIF.freeparametercoordinate));
                hatx.y=soln.solverinfo.parameters(sol.solverinfo.equilibriumcoord);
                hatx.x=0;
                hatx.arcarg=OCMATCONT.HE.arcarg(end);
                hatx.linearization=jacobian(dgObj0tmp,hatx);
                ocEP=dynprimitive(hatx);
                ocResultStruct.NonadmissibleSolution=ocasymptotic(octrajectory(soln,dgObj0tmp),ocEP);
            elseif any(strcmp(fieldname,{'indifferencesolution4ep_ft','indifferencesolutionep'}))
                parnew=soln.modelparameter;
                parnew(OCMATINDIF.parameterindex)=soln.solverinfo.parameters(OCMATINDIF.parametervaluecoord);
                dgObj0tmp=changeparametervalue(dgObj0,parnew);
                J=cell(1,OCMATINDIF.indifferenceorder);
                saddlepoint=cell(1,OCMATINDIF.indifferenceorder);
                for ii=1:OCMATINDIF.indifferenceorder
                    if OCMATINDIF.ocasymptotic(ii)
                        equcoord=soln.solverinfo.equilibriumcoord{ii};
                    else
                        equcoord=[];
                    end
                    if ~isempty(equcoord)
                        hatx.y=soln.solverinfo.parameters(equcoord);
                        hatx.x=0;
                        hatx.arcarg=OCMATCONT.HE.arcarg(OCMATINDIF.arccoord{ii}(end));
                        J{ii}=jacobian(dgObj0tmp,hatx);
                        saddlepoint{ii}=hatx.y;
                    end
                end
                ocResultStruct.NonadmissibleSolution=sol2ocmultipath(soln,OCMATINDIF,saddlepoint,J);
                parnew=sol.modelparameter;
                parnew(OCMATINDIF.parameterindex)=sol.solverinfo.parameters(OCMATINDIF.parametervaluecoord);
                dgObj0tmp=changeparametervalue(dgObj0,parnew);
                J=cell(1,OCMATINDIF.indifferenceorder);
                saddlepoint=cell(1,OCMATINDIF.indifferenceorder);
                for ii=1:OCMATINDIF.indifferenceorder
                    if OCMATINDIF.ocasymptotic(ii)
                        equcoord=sol.solverinfo.equilibriumcoord{ii};
                    else
                        equcoord=[];
                    end
                    if ~isempty(equcoord)
                        hatx.y=sol.solverinfo.parameters(equcoord);
                        hatx.x=0;
                        hatx.arcarg=OCMATCONT.HE.arcarg(OCMATINDIF.arccoord{ii}(end));
                        J{ii}=jacobian(dgObj0tmp,hatx);
                        saddlepoint{ii}=hatx.y;
                    end
                end
                ocResultStruct.ExtremalSolution=sol2ocmultipath(sol,OCMATINDIF,saddlepoint,J);
            end
        else
            if ~isempty(soltarget)
                if any(strcmp(fieldname,{'indifferencesolution','indifferencedistribution','indifferencesolution4ae2ftae'}))
                    switch soltarget.solver
                        case 'gbvp4c'
                            ocResultStruct.ExtremalSolution=sol2ocmultipath(soltarget,OCMATINDIF,OCMATINDIF.limitset,OCMATINDIF.linearization);
                        otherwise
                            ocResultStruct.ExtremalSolution=sol2ocmultipath(soltarget,OCMATINDIF,OCMATINDIF.saddlepoint,OCMATINDIF.linearization);
                    end
                elseif strcmp(fieldname,'indifferencesolution4per')
                    ocResultStruct.ExtremalSolution=sol2ocmultipath(soltarget,OCMATINDIF,OCMATINDIF.limitset,OCMATINDIF.linearization);
                elseif strcmp(fieldname,'indifferencesolution4emf')
                    ocResultStruct.ExtremalSolution=sol2ocmultipath(soltarget,OCMATINDIF,soltarget.solverinfo.saddlepoint,OCMATINDIF.linearization);
                elseif strcmp(fieldname,'indifferencesolution4ft')
                    ocResultStruct.ExtremalSolution=sol2ocmultipath(soltarget);
                elseif strcmp(fieldname,'indifferencesolutionep')
                    dgObj0tmp=changeparametervalue(dgObj0,OCMATINDIF.freeparameterindex,soltarget.solverinfo.parameters(OCMATINDIF.freeparametercoordinate));
                    hatx.y=soltarget.solverinfo.parameters(soltarget.solverinfo.equilibriumcoord);
                    hatx.x=0;
                    hatx.arcarg=OCMATCONT.HE.arcarg(end);
                    hatx.linearization=jacobian(dgObj0tmp,hatx);
                    ocEP=dynprimitive(hatx);
                    ocResultStruct.ExtremalSolution=ocasymptotic(octrajectory(soltarget,dgObj0tmp),ocEP);
                elseif strcmp(fieldname,'indifferencesolution4ep_ft')
                    parnew=soltarget.modelparameter;
                    parnew(OCMATINDIF.parameterindex)=soltarget.solverinfo.parameters(OCMATINDIF.parametervaluecoord);
                    dgObj0tmp=changeparametervalue(dgObj0,parnew);
                    J=cell(1,OCMATINDIF.indifferenceorder);
                    saddlepoint=cell(1,OCMATINDIF.indifferenceorder);
                    for ii=1:OCMATINDIF.indifferenceorder
                        if OCMATINDIF.ocasymptotic(ii)
                            equcoord=soltarget.solverinfo.equilibriumcoord{ii};
                        else
                            equcoord=[];
                        end
                        if ~isempty(equcoord)
                            hatx.y=soltarget.solverinfo.parameters(equcoord);
                            hatx.x=0;
                            hatx.arcarg=OCMATCONT.HE.arcarg(OCMATINDIF.arccoord{ii}(end));
                            J{ii}=jacobian(dgObj0tmp,hatx);
                            saddlepoint{ii}=hatx.y;
                        end
                    end
                    ocResultStruct.ExtremalSolution=sol2ocmultipath(soltarget,OCMATINDIF,saddlepoint,J);
                end
            else
                if any(strcmp(fieldname,{'indifferencesolution','indifferencedistribution','indifferencesolution4ae2ftae'}))
                    switch ocExStruct(end).solver
                        case 'gbvp4c'
                            ocResultStruct.ExtremalSolution=sol2ocmultipath(ocExStruct(end),OCMATINDIF,OCMATINDIF.limitset,OCMATINDIF.linearization);
                        otherwise
                            ocResultStruct.ExtremalSolution=sol2ocmultipath(ocExStruct(end),OCMATINDIF,OCMATINDIF.limitset,OCMATINDIF.linearization);
                    end
                elseif strcmp(fieldname,'indifferencesolution4per')
                    ocResultStruct.ExtremalSolution=sol2ocmultipath(ocExStruct(end),OCMATINDIF,OCMATINDIF.limitset,OCMATINDIF.linearization);
                elseif strcmp(fieldname,'indifferencesolution4emf')
                    ocResultStruct.ExtremalSolution=sol2ocmultipath(ocExStruct(end),OCMATINDIF,ocExStruct(end).solverinfo.saddlepoint,OCMATINDIF.linearization);
                elseif strcmp(fieldname,'indifferencesolution4ft')
                    ocResultStruct.ExtremalSolution=sol2ocmultipath(ocExStruct(end));
                elseif strcmp(fieldname,'indifferencesolutionep')
                    dgObj0tmp=changeparametervalue(dgObj0,OCMATINDIF.freeparameterindex,ocExStruct(end).solverinfo.parameters(OCMATINDIF.freeparametercoordinate));
                    hatx.y=ocExStruct(end).solverinfo.parameters(ocExStruct(end).solverinfo.equilibriumcoord);
                    hatx.x=0;
                    hatx.arcarg=OCMATCONT.HE.arcarg(end);
                    hatx.linearization=jacobian(dgObj0tmp,hatx);
                    ocEP=dynprimitive(hatx);
                    ocResultStruct.ExtremalSolution=ocasymptotic(octrajectory(ocExStruct(end),dgObj0tmp),ocEP);
                elseif strcmp(fieldname,'indifferencesolution4ep_ft')
                    parnew=ocExStruct(end).modelparameter;
                    parnew(OCMATINDIF.parameterindex)=ocExStruct(end).solverinfo.parameters(OCMATINDIF.parametervaluecoord);
                    dgObj0tmp=changeparametervalue(dgObj0,parnew);
                    J=cell(1,OCMATINDIF.indifferenceorder);
                    saddlepoint=cell(1,OCMATINDIF.indifferenceorder);
                    for ii=1:OCMATINDIF.indifferenceorder
                        if OCMATINDIF.ocasymptotic(ii)
                            equcoord=ocExStruct(end).solverinfo.equilibriumcoord{ii};
                        else
                            equcoord=[];
                        end
                        if ~isempty(equcoord)
                            hatx.y=ocExStruct(end).solverinfo.parameters(equcoord);
                            hatx.x=0;
                            hatx.arcarg=OCMATCONT.HE.arcarg(OCMATINDIF.arccoord{ii}(end));
                            J{ii}=jacobian(dgObj0tmp,hatx);
                            saddlepoint{ii}=hatx.y;
                        end
                    end
                    ocResultStruct.ExtremalSolution=sol2ocmultipath(ocExStruct(end),OCMATINDIF,saddlepoint,J);
                end
                if strcmp(fieldname,'indifferencesolution4ft')
                    ocResultStruct.NonadmissibleSolution=octrajectory([]);
                else
                    ocResultStruct.NonadmissibleSolution=ocasymptotic([]);
                end
            end
        end
        if any(strcmp(fieldname,{'indifferencesolution','indifferencesolution4emf','indifferencesolution4per','indifferencedistribution','indifferencesolution4ft','indifferencesolution4ae2ftae','indifferencesolution4ep_ft','indifferencesolutionep'}))
            ocResultStruct.IndifferencePointCurve=occurve(indifferencepointStruct);
        end
        ocResultStruct.ContinuationSolution=ocExStruct;
        ocResultStruct.ContinuationClassification=fieldname;
        ocResultStruct.ContinuationInformation=sout;
        if any(strcmp(fieldname,{'indifferencesolution4ft'})) && OCMATINDIF.continuationtype==3
            dgObj0tmp=changeparametervalue(dgObj0,OCMATINDIF.parameterindex,ocExStruct(end).solverinfo.parameters(OCMATINDIF.parametercoord));
            ocResultStruct.LastModel=dgObj0tmp;
        elseif any(strcmp(fieldname,{'indifferencesolution4ep_ft'}))
            ocResultStruct.LastModel=dgObj0tmp;
        end
        %ocResultStruct.LimitSet=ocEP;
        ocResultFieldName='Continuation';

        resetglobal(fieldname,OCMATCONTORIGINAL,OCMATINDIFORIGINAL,OCBVPORIGINAL);
        
    case {'modelequilibrium'}
        global OCMATFINITCONT

        % store value of the global variables OCMATCONT and OCMATAE to
        % temporary variables, in an intermediate step the global variables
        % are changed to values stored in the continuation files and reset
        % afterwards
        if ~isempty(OCMATFINITCONT)
            OCMATFINITCONTORIGINAL=OCMATFINITCONT;
        else
            OCMATFINITCONTORIGINAL=[];
        end
        % the file names containing continuation information are build by
        % the model specific basic names and the continuation
        % classification, e.g. 4extremal2ep
        resultfile=[OCMATFINITCONT.basicresultfilename '4modelequilibrium.mat'];
        globalvarfile=[OCMATFINITCONT.basicglobalvarfilename '4modelequilibrium.mat'];

        % test if files exist
        if ~exist(resultfile,'file')
            ocmatmsg('Continuation file ''%s'' does not exist\n',resultfile)
            return
        end
        if ~exist(globalvarfile,'file')
            ocmatmsg('Continuation file ''%s'' does not exist\n',globalvarfile)
            return
        end

        % test if data file is older than 12 hours in that case ask if user
        % want to proceed or interrupt
        info=dir(resultfile);
        if abs(info.datenum-datenum(now))>0.5
            ocmatmsg('Data file was already generated at %s.\n',datestr(info.datenum))
            answer=input('Interrupt storing?  y/(n): ','s');
            if isempty(answer)
                % default value 'y'
                answer='n';
            end
            if strcmpi(answer,'y')
                disp('Return without storing.')
                return
            end
        end

        % load data of information files into workspace
        load(globalvarfile)
        load(resultfile)
        [numrow numcol]=size(xout);
        switch fieldname
            case 'modelequilibrium'
                ocC.y=xout(1:MODELINFO.MATCONTEQUILIBRIUM.nphase,1:numcol);
                ocC.arcarg=OCMATFINITCONT.arcarg(ones(1,numcol));
                ocC.arcposition=[1 numcol];
                ocC.modelname=modelname(dgObj);
                ocC.modelparameter=MODELINFO.MATCONTEQUILIBRIUM.P0;
                ocC.userinfo.varyparameterindex=MODELINFO.MATCONTEQUILIBRIUM.ActiveParams;
                ocC.userinfo.varyparametervalue=xout(MODELINFO.MATCONTEQUILIBRIUM.nphase+1,1:numcol);
        end
        ocC.userinfo.tangent=vout;
        ocResultStruct.ContinuationInformation=sout;
        ocResultStruct.ContinuationSolution=occurve(ocC);

        ocResultStruct.ContinuationClassification=fieldname;
        ocResultFieldName='MatContContinuation';

    otherwise
        ocResultStruct=fieldvalue;
        ocResultFieldName=fieldname;
end

function [globalvarfile,resultfile,varargout]=loadglobal(dgObj,fieldname)
switch fieldname
    case {'extremaldg2ep','extremal2ep'}
        global OCMATCONT OCMATAE OCBVP
        % store value of the global variables OCMATCONT and OCMATAE to
        % temporary variables, in an intermediate step the global variables
        % are changed to values stored in the continuation files and reset
        % afterwards
        if ~isempty(OCMATCONT)
            OCMATCONTORIGINAL=OCMATCONT;
        else
            OCMATCONTORIGINAL=[];
        end
        if ~isempty(OCMATAE)
            OCMATAEORIGINAL=OCMATAE;
        else
            OCMATAEORIGINAL=[];
        end
        if ~isempty(OCBVP)
            OCBVPORIGINAL=OCBVP;
        else
            OCBVPORIGINAL=[];
        end

        % the file names containing continuation information are build by
        % the model specific basic names and the continuation
        % classification, e.g. 4extremal2ep
        if ~isempty(OCMATAE)
            resultfile=[OCMATAE.basicresultfilename '4' fieldname '.mat'];
            globalvarfile=[OCMATAE.basicglobalvarfilename '4' fieldname '.mat'];
        end

        if isempty(OCMATAE) || ~exist(resultfile,'file')
            resultfile=fullfile(getocmatpath,getocmatfolder(dgObj,'userdata'),['SaveIntermediateResults4' fieldname '.mat']);
        end
        if isempty(OCMATAE) || ~exist(globalvarfile,'file')
            globalvarfile=fullfile(getocmatpath,getocmatfolder(dgObj,'userdata'),['SaveIntermediateResultsGlobalVariable4' fieldname '.mat']);
        end
        % test if files exist
        if ~exist(resultfile,'file')
            ocmatmsg('Continuation file ''%s'' does not exist\n',resultfile)
            return
        end
        if ~exist(globalvarfile,'file')
            ocmatmsg('Continuation file ''%s'' does not exist\n',globalvarfile)
            return
        end

        % test if data file is older than 12 hours in that case ask if user
        % want to proceed or interrupt
        info=dir(resultfile);
        if abs(info.datenum-datenum(now))>0.5
            ocmatmsg('Data file was already generated at %s.\n',datestr(info.datenum))
            answer=input('Interrupt storing?  y/(n): ','s');
            if isempty(answer)
                % default value 'y'
                answer='n';
            end
            if strcmpi(answer,'y')
                disp('Return without storing.')
                return
            end
        end
        % load data of information files into workspace
        load(globalvarfile)
        load(resultfile)

        % change actual global variables by the values stored in the
        % continuation file as MODELINFO.OCMATCONT and MODELINFO.OCMATAE.
        % If the storing process is interrupted before the global variables
        % are reset to its original value the user is informed that s/he
        % has to newly initialize an actual continuation process.
        ocmatmsg('Global variable ''OCMATCONT'' changed.\n')
        OCMATCONT=MODELINFO.OCMATCONT;
        OCMATCONT.meshadaptationflag=1;
        ocmatmsg('Global variable ''OCMATAE'' changed.\n')
        OCMATAE=MODELINFO.OCMATAE;
        ocmatmsg('Global variable ''OCBVP'' changed.\n')
        OCBVP=MODELINFO.OCBVP;
        varargout{1}=OCMATCONTORIGINAL;
        varargout{2}=OCMATAEORIGINAL;
        varargout{3}=OCBVPORIGINAL;

    case {'indifferencesolution'}
        global OCMATCONT OCBVP OCMATINDIF
        % store value of the global variables OCMATCONT and OCMATAE to
        % temporary variables, in an intermediate step the global variables
        % are changed to values stored in the continuation files and reset
        % afterwards
        if ~isempty(OCMATCONT)
            OCMATCONTORIGINAL=OCMATCONT;
        else
            OCMATCONTORIGINAL=[];
        end
        if ~isempty(OCBVP)
            OCBVPORIGINAL=OCBVP;
        else
            OCBVPORIGINAL=[];
        end
        if ~isempty(OCMATINDIF)
            OCMATINDIFORIGINAL=OCMATINDIF;
        else
            OCMATINDIFORIGINAL=[];
        end

        % the file names containing continuation information are build by
        % the model specific basic names and the continuation
        % classification, e.g. 4extremal2ep
        if ~isempty(OCMATINDIF)
            resultfile=[OCMATINDIF.basicresultfilename '4' fieldname '.mat'];
            globalvarfile=[OCMATINDIF.basicglobalvarfilename '4' fieldname '.mat'];
        else
            resultfile=fullfile(getocmatpath,getocmatfolder(dgObj0,'userdata'),['SaveIntermediateResults4' fieldname '.mat']);
            globalvarfile=fullfile(getocmatpath,getocmatfolder(dgObj0,'userdata'),['SaveIntermediateResultsGlobalVariable4' fieldname '.mat']);
        end

        % test if files exist
        if ~exist(resultfile,'file')
            ocmatmsg('Continuation file ''%s'' does not exist\n',resultfile)
            return
        end
        if ~exist(globalvarfile,'file')
            ocmatmsg('Continuation file ''%s'' does not exist\n',globalvarfile)
            return
        end

        % test if data file is older than 12 hours in that case ask if user
        % want to proceed or interrupt
        info=dir(resultfile);
        if abs(info.datenum-datenum(now))>0.5
            ocmatmsg('Data file has been generated at %s.\n',datestr(info.datenum))
            answer=input('Interrupt storing?  y/(n): ','s');
            if isempty(answer)
                % default value 'y'
                answer='n';
            end
            if strcmpi(answer,'y')
                disp('Return without storing.')
                return
            end
        end
        % load data of information files into workspace
        load(globalvarfile)
        load(resultfile)

        % change actual global variables by the values stored in the
        % continuation file as MODELINFO.OCMATCONT and MODELINFO.OCMATAE.
        % If the storing process is interrupted before the global variables
        % are reset to its original value the user is informed that s/he
        % has to newly initialize an actual continuation process.
        ocmatmsg('Global variable ''OCMATCONT'' changed.\n')
        OCMATCONT=MODELINFO.OCMATCONT;
        OCMATCONT.meshadaptationflag=1;
        ocmatmsg('Global variable ''OCMATINDIF'' changed.\n')
        OCMATINDIF=MODELINFO.OCMATINDIF;
        ocmatmsg('Global variable ''OCBVP'' changed.\n')
        OCBVP=MODELINFO.OCBVP;
        varargout{1}=OCMATCONTORIGINAL;
        varargout{2}=OCMATINDIFORIGINAL;
        varargout{3}=OCBVPORIGINAL;
        
end

function resetglobal(fieldname,varargin)
switch fieldname
    case {'extremaldg2ep','extremal2ep'}
        OCMATCONTORIGINAL=varargin{1};
        OCMATAEORIGINAL=varargin{2};
        OCBVPORIGINAL=varargin{3};
        % reset global variables to its original value or delete it if no
        % global variables existed
        if ~isempty(OCMATCONTORIGINAL)
            ocmatmsg('Global variable ''OCMATCONT'' reset.\n')
            OCMATCONT=OCMATCONTORIGINAL;
        else
            clear global OCMATCONT
        end
        if ~isempty(OCMATAEORIGINAL)
            ocmatmsg('Global variable ''OCMATAE'' reset.\n')
            OCMATAE=OCMATAEORIGINAL;
        else
            clear global OCMATAE
        end
        if ~isempty(OCBVPORIGINAL)
            ocmatmsg('Global variable ''OCBVP'' reset.\n')
            OCBVP=OCBVPORIGINAL;
        else
            clear global OCBVP
        end

    case {'indifferencesolution'}
        OCMATCONTORIGINAL=varargin{1};
        OCMATINDIFORIGINAL=varargin{2};
        OCBVPORIGINAL=varargin{3};
        % reset global variables to its original value or delete it if no
        % global variables existed
        if ~isempty(OCMATCONTORIGINAL)
            ocmatmsg('Global variable ''OCMATCONT'' reset.\n')
            OCMATCONT=OCMATCONTORIGINAL;
        else
            clear global OCMATCONT
        end
        if ~isempty(OCMATINDIFORIGINAL)
            ocmatmsg('Global variable ''OCMATINDIF'' reset.\n')
            OCMATINDIF=OCMATINDIFORIGINAL;
        else
            clear global OCMATINDIF
        end
        if ~isempty(OCBVPORIGINAL)
            ocmatmsg('Global variable ''OCBVP'' reset.\n')
            OCBVP=OCBVPORIGINAL;
        else
            clear global OCBVP
        end
end

function [solt,sole,soln,sol]=findspecificsolution(sout)
% write violation information into violationinfo field of soln
solt=[]; % solution at target point
sole=[]; % last solution
sol=[];
soln=[];
if isempty(sout)
    return
end
counter=numel(sout);
sole=sout(counter).data.sol;
while counter>2
    counter=counter-1;
    if strcmp(sout(counter).label,'NAS')
        sol=sout(counter).data.sol;
        soln=sout(counter).data.soln;
        sout(counter).data=rmfield(sout(counter).data,{'sol','soln'});
        soln.violationinfo=sout(counter).data.infoS;
    end
    if strcmp(sout(counter).label,'HTV')
        solt=sout(counter).data.sol;
    end
end
