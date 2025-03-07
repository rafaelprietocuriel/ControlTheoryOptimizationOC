function odeObj=store(odeObj,varargin)
%
% STORE results
%
% ODEOBJ=STORE(ODEOBJ,'CONTINUATION') store results of the last
% continuation process in the results of the oc model ODEOBJ. The identifier
% 'CONTINUATION' can also be omitted.
%
% ODEOBJ=STORE(ODEOBJ,'DATAFILE.MAT') store data given by the MAT file DATAFILE.MAT as
% a result of a continuation process.
%
% ODEOBJ=STORE(ODEOBJ,OCELEMENT) stores the object OCELEMENT in the results of the oc
% object ODEOBJ (OCELEMENT can be a cell aray of objects). Possible classes for
% OCELEMENT are
%        dynprimitive: equilibrium or limit cycle
%        ocasymptotic: asymptotic solutions
%
% The default fieldnames of the ocResults are
%       ocResults.Equilibrium: for a dynprimitive being an equilibrium
%       ocResults.LimitCycle: for a dynprimitive being a limit cycle
%       ocResults.ExtremalSolution: for an ocasymptotic
%       ocResults.NonAdmissibleSolution: for an ocasymptotic violating the
%                                        models constraints
%
% ODEOBJ=STORE(...,'REPLACE',VAL) forces to replace an already calculated
% object by the new one if VAL='on' (default it is set 'off'). If the
% object does not already exist it is appended.
%
% ODEOBJ=STORE(...,'APPEND',VAL) forces to append the object to an already
% existing cell array of objects or creates a new cell array.
%
% ODEOBJ=STORE(...,'FIELDNAME',VAL) forces to store the object at the
% ocResults under the fieldname VAL, if the object is a cell array, VAL is
% either a string or a cell array of strings of the same length as VAL.
%
% [ODEOBJ,OCELEMENT]=STORE(ODEOBJ,...) returns the stored elements OCELEMENT.
%
% [ODEOBJ,OCELEMENT,FIELDNAMES]=STORE(ODEOBJ,...) returns corresponding
% fieldnames.


if isempty(odeObj) || nargin==1
    return
end

fieldname=[];
fieldvalue=[];
if ~ischar(varargin{1}) && nargin==2
    if isdynprimitive(varargin{1}) ||  isocasymptotic(varargin{1}) || isoctrajectory(varargin{1}) || isoccurve(varargin{1})
        fieldname=class(varargin{1});
    elseif iscell(varargin{1})
        for ii=1:numel(varargin{1})
            store(odeObj,varargin{1}{ii},varargin{2:end});
        end

        if ~nargout
            assignin('caller',inputname(1),odeObj);
        else
            varargout{1}=odeObj;
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
[ocStruct,ocResultFieldName]=generateelement(odeObj,fieldname,fieldvalue);
if isempty(ocResultFieldName)
    return
end
if isfield(odeObj.Result,ocResultFieldName)
    odeObj.Result.(ocResultFieldName){end+1}=ocStruct;
else
    odeObj.Result(1).(ocResultFieldName){1}=ocStruct;
end

if ~nargout
    assignin('caller',inputname(1),odeObj);
else
    varargout{1}=odeObj;
end

%%
function [ocResultStruct,ocResultFieldName]=generateelement(odeObj,fieldname,fieldvalue)
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
    case {'saddlepath2ep'}
        global OCMATCONT OCMATAE OCBVP

        [globalvarfile,resultfile,OCMATCONTORIGINAL,OCMATAEORIGINAL,OCBVPORIGINAL]=loadglobal(odeObj,fieldname);

        % load data of information files into workspace
        load(globalvarfile)
        load(resultfile)

        mainfunch=feval(str2func(fieldname));

        % determine if non admissible solution occurs
        [soltarget sollast soln sol]=findspecificsolution(sout);
        counter=0;
        for ii=1:numel(bvpout)
            try
                counter=counter+1;
                DataAdaptation(bvpout(ii).tmesh);
                ocExStruct(counter)=mainfunch{22}(bvpout(ii).tmesh,bvpout(ii).coeff,bvpout(ii).tangent);
                slMfStruct.y(OCMATCONT.DOMAINDDATA(OCMATCONT.HE.arcindex(1)).eqcoord,counter)=ocExStruct(ii).y(OCMATCONT.DOMAINDDATA(OCMATCONT.HE.arcindex(1)).eqcoord,1);
                slMfStruct.arcarg(counter)=OCMATCONT.HE.arcarg(1);
                slMfStruct.userinfo.continuationparameter(counter)=bvpout(ii).coeff(end);
            end
        end
        slMfStruct.arcposition=[];
        ocEP=dynprimitive(OCMATAE.saddlepoint,OCMATCONT.HE.arcarg(end),OCMATAE.linearization);
        if ~isempty(soln)
            ocResultStruct.NonadmissibleSolution=ocasymptotic(octrajectory(soln,odeObj),ocEP);
            ocResultStruct.ExtremalSolution=ocasymptotic(octrajectory(sol,odeObj),ocEP);
            matlabsol=transform2nativematlab(ocExStruct(end).x,ocExStruct(end).solverinfo.coeff,ocExStruct(end).modelparameter);
            if ~isempty(matlabsol)
                ocResultStruct.ExtremalSolutionMatlabFormat=matlabsol;
            end
        else
            if ~isempty(soltarget)
                ocResultStruct.ExtremalSolution=ocasymptotic(octrajectory(soltarget,odeObj),ocEP);
                matlabsol=transform2nativematlab(soltarget.solverinfo.tmesh,soltarget.solverinfo.coeff,soltarget.modelparameter);
                if ~isempty(matlabsol)
                    ocResultStruct.ExtremalSolutionMatlabFormat=matlabsol;
                end
            else
                DataAdaptation(ocExStruct(end).solverinfo.tmesh);
                ocResultStruct.ExtremalSolution=ocasymptotic(octrajectory(ocExStruct(end),odeObj),ocEP);
                matlabsol=transform2nativematlab(ocExStruct(end).solverinfo.tmesh,ocExStruct(end).solverinfo.coeff,ocExStruct(end).modelparameter);
                if ~isempty(matlabsol)
                    ocResultStruct.ExtremalSolutionMatlabFormat=matlabsol;
                end
            end
            ocResultStruct.NonadmissibleSolution=ocasymptotic([]);
        end
        if strcmp(fieldname,'extremal2ep')
            ocResultStruct.SliceManifold=occurve(slMfStruct);
        end
        ocResultStruct.ContinuationSolution=ocExStruct;
        ocResultStruct.ContinuationClassification=fieldname;
        ocResultStruct.ContinuationInformation=sout;
        ocResultStruct.LimitSet=ocEP;
        ocResultFieldName='Continuation';

        resetglobal(fieldname,OCMATCONTORIGINAL,OCMATAEORIGINAL,OCBVPORIGINAL);
        
    case {'odesolution','userbvp','odesolution4t','odesolution4s','odesolution4p'}
        global OCMATCONT OCMATODEBVP OCBVP
        [globalvarfile,resultfile,OCMATCONTORIGINAL,OCMATODEBVPORIGINAL,OCBVPORIGINAL]=loadglobal(odeObj,fieldname);

        % load data of information files into workspace
        load(globalvarfile)
        load(resultfile)

        mainfunch=feval(str2func(fieldname));
        % determine if non admissible solution occurs
        [soltarget sollast soln sol]=findspecificsolution(sout);

        counter=1;
        for ii=1:numel(bvpout)
            try
                DataAdaptation(bvpout(ii).tmesh);
                ocExStruct(counter)=mainfunch{22}(bvpout(ii).tmesh,bvpout(ii).coeff,bvpout(ii).tangent);
                counter=counter+1;
            end
        end

        if ~isempty(soltarget)
            sol=soltarget;
            if any(strcmp(fieldname,{'odesolution4p','odesolution4s'}))
                if ~isempty(OCMATODEBVP.changeparameterindex)
                    odeObjtmp=changeparametervalue(odeObj,OCMATODEBVP.changeparameterindex,sol.solverinfo.parameters(sol.solverinfo.changeparametercoord));
                else
                    odeObjtmp=odeObj;
                end
            else
                odeObjtmp=odeObj;
            end
            ocResultStruct.LastSolution=octrajectory(sol,odeObjtmp);
            ocResultStruct.NonadmissibleSolution=octrajectory([]);
        end

        if ~isempty(sollast)
            sol=sollast;
            if any(strcmp(fieldname,{'odesolution4p','odesolution4s'}))
                if ~isempty(OCMATODEBVP.changeparameterindex)
                    odeObjtmp=changeparametervalue(odeObj,OCMATODEBVP.changeparameterindex,sol.solverinfo.parameters(sol.solverinfo.changeparametercoord));
                else
                    odeObjtmp=odeObj;
                end
            else
                odeObjtmp=odeObj;
            end
            ocResultStruct.NonadmissibleSolution=octrajectory(sol,odeObjtmp);
            sol=ocExStruct(end);
            if any(strcmp(fieldname,{'odesolution4p','odesolution4s'}))
                if ~isempty(OCMATODEBVP.changeparameterindex)
                    odeObjtmp=changeparametervalue(odeObj,OCMATODEBVP.changeparameterindex,sol.solverinfo.parameters(sol.solverinfo.changeparametercoord));
                else
                    odeObjtmp=odeObj;
                end
            else
                odeObjtmp=odeObj;
            end
            ocResultStruct.LastSolution=octrajectory(sol,odeObjtmp);
        end
        if ~isempty(soltarget)
            sol=soltarget;
            if any(strcmp(fieldname,{'odesolution4p','odesolution4s'}))
                if ~isempty(OCMATODEBVP.changeparameterindex)
                    odeObjtmp=changeparametervalue(odeObj,OCMATODEBVP.changeparameterindex,sol.solverinfo.parameters(sol.solverinfo.changeparametercoord));
                else
                    odeObjtmp=odeObj;
                end
            else
                odeObjtmp=odeObj;
            end
            ocResultStruct.LastSolution=octrajectory(sol,odeObjtmp);
            ocResultStruct.NonadmissibleSolution=octrajectory([]);
        end
        %odeObjtmp=changeparametervalue(odeObj,OCMATODEBVP.varyparameterindex,sol.solverinfo.continuationparameter);
        %ocResultStruct.ContinuationParameter=OCMATODEBVP.varyparameterindex;
        %ocResultStruct.LastSolution=octrajectory(ocExStruct(end));
        ocResultStruct.ContinuationSolution=ocExStruct;
        ocResultStruct.ContinuationClassification=fieldname;
        ocResultStruct.ContinuationInformation=sout;
        if any(strcmp(fieldname,{'odesolution4p','odesolution4s'}))
            if ~isempty(OCMATODEBVP.changeparameterindex)
                ocResultStruct.LastModel=odeObjtmp;
            end
        end
        %
        ocResultFieldName='Continuation';
        
        resetglobal(fieldname,OCMATCONTORIGINAL,OCMATODEBVPORIGINAL,OCBVPORIGINAL);
        
        
    case {'modelequilibrium','modelequilibriummf'}
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
                ocC.modelname=modelname(odeObj);
                ocC.modelparameter=MODELINFO.MATCONTEQUILIBRIUM.P0;
                ocC.userinfo.varyparameterindex=MODELINFO.MATCONTEQUILIBRIUM.ActiveParams;
                ocC.userinfo.varyparametervalue=xout(MODELINFO.MATCONTEQUILIBRIUM.nphase+1,1:numcol);
            case 'modelequilibriummf'
                ocC.y=xout(1:MODELINFO.MATCONTEQUILIBRIUM.nphase+1,1:numcol);
                ocC.arcarg=OCMATFINITCONT.arcarg(ones(1,numcol));
                ocC.arcposition=[1 numcol];
                ocC.modelname=modelname(odeObj);
                ocC.modelparameter=MODELINFO.MATCONTEQUILIBRIUM.P0;
                ocC.userinfo.varyparameterindex=MODELINFO.MATCONTEQUILIBRIUM.ActiveParams;
                ocC.userinfo.varyparametervalue=xout(MODELINFO.MATCONTEQUILIBRIUM.nphase+1,1:numcol);
        end
        ocC.userinfo.tangent=vout;
        ocResultStruct.ContinuationInformation=sout;
        ocResultStruct.ContinuationSolution=occurve(ocC);

        ocResultStruct.ContinuationClassification=fieldname;
        ocResultFieldName='MatContContinuation';


    case {'limitcycle','limitcycleuser','limitcycleuserI'}
        global OCMATCONT OCMATLC OCBVP

        [globalvarfile,resultfile,OCMATCONTORIGINAL,OCMATAEORIGINAL,OCBVPORIGINAL]=loadglobal(odeObj,fieldname);

        % load data of information files into workspace
        load(globalvarfile)
        load(resultfile)

        mainfunch=feval(str2func(fieldname));

        % determine if non admissible solution occurs
        [soltarget sollast soln sol]=findspecificsolution(sout);
        counter=0;
        OCBVP.F=[];
        for ii=1:numel(bvpout)
            try
                counter=counter+1;
                DataAdaptation(bvpout(ii).tmesh);
                ocExStruct(counter)=mainfunch{22}(bvpout(ii).tmesh,bvpout(ii).coeff,bvpout(ii).tangent);
            catch
                lasterr
            end
        end
        if ~isempty(soln)
            if any(strcmp(fieldname,{'limitcycle','limitcycleuser'}))
                odeObjtmp=changeparametervalue(odeObj,OCMATLC.varyparameterindex,soln.octrajectory.solverinfo.continuationparameter);
            else
                odeObjtmp=odeObj;
            end
            ocResultStruct.NonadmissibleLimitCycle=dynprimitive(soln,odeObjtmp);
            ocResultStruct.LimitCycle=dynprimitive(sol,odeObjtmp);
            trajectory=ocExStruct(end).octrajectory;
            matlabsol=transform2nativematlab(trajectory.solverinfo.tmesh,trajectory.solverinfo.coeff,trajectory.modelparameter);
            if ~isempty(matlabsol)
                ocResultStruct.LimitCycleMatlabFormat=matlabsol;
            end
        else
            if ~isempty(soltarget)
                if any(strcmp(fieldname,{'limitcycle','limitcycleuser'}))
                    odeObjtmp=changeparametervalue(odeObj,OCMATLC.varyparameterindex,soltarget.octrajectory.solverinfo.continuationparameter);
                else
                    odeObjtmp=odeObj;
                end
                ocResultStruct.LimitCycle=dynprimitive(soltarget,odeObjtmp);
                trajectory=soltarget.octrajectory;
                matlabsol=transform2nativematlab(trajectory.solverinfo.tmesh,trajectory.solverinfo.coeff,trajectory.modelparameter);
                if ~isempty(matlabsol)
                    ocResultStruct.LimitCycleMatlabFormat=matlabsol;
                end
            else
                DataAdaptation(ocExStruct(end).octrajectory.solverinfo.tmesh);
                if any(strcmp(fieldname,{'limitcycle','limitcycleuser'}))
                    odeObjtmp=changeparametervalue(odeObj,OCMATLC.varyparameterindex,ocExStruct(end).octrajectory.solverinfo.continuationparameter);
                else
                    odeObjtmp=odeObj;
                end
                ocResultStruct.LimitCycle=dynprimitive(sollast);
                trajectory=ocExStruct(end).octrajectory;
                matlabsol=transform2nativematlab(trajectory.solverinfo.tmesh,trajectory.solverinfo.coeff,trajectory.modelparameter);
                if ~isempty(matlabsol)
                    ocResultStruct.LimitCycleMatlabFormat=matlabsol;
                end
            end
            ocResultStruct.NonadmissibleLimitCycle=dynprimitive([]);
        end
        ocResultStruct.ContinuationSolution=ocExStruct;
        ocResultStruct.ContinuationClassification=fieldname;
        ocResultStruct.ContinuationInformation=sout;
        ocResultStruct.LastModel=odeObjtmp;
        ocResultFieldName='Continuation';

        resetglobal(fieldname,OCMATCONTORIGINAL,OCMATAEORIGINAL,OCBVPORIGINAL);

    case 'modellimitcycle'
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
        resultfile=[OCMATFINITCONT.basicresultfilename '4modellimitcycle.mat'];
        globalvarfile=[OCMATFINITCONT.basicglobalvarfilename '4modellimitcycle.mat'];

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
        totalmesh=fout(1:MODELINFO.MATCONTLIMITCYCLE.ntst+1,:);
        if isfield(MODELINFO.MATCONTLIMITCYCLE,'multipliers') && ~isempty(MODELINFO.MATCONTLIMITCYCLE.multipliers)
            linearization=fout(MODELINFO.MATCONTLIMITCYCLE.ntst+MODELINFO.MATCONTLIMITCYCLE.nphase+2:end,:);
        else
            linearization=[];
        end
        ntstcolp1=MODELINFO.MATCONTLIMITCYCLE.ntstcol+1;
        nonemptyidx=find(sum(abs(xout)));
        modelnme=modelname(odeObj);
        counter=0;
        for ii=nonemptyidx
            counter=counter+1;
            trj(counter).y=reshape(xout(MODELINFO.MATCONTLIMITCYCLE.coords,ii),MODELINFO.MATCONTLIMITCYCLE.nphase,MODELINFO.MATCONTLIMITCYCLE.tps);
            msh=totalmesh(:,ii).';
            cumsum([msh(1:end-1);repmat(diff(msh)/MODELINFO.MATCONTLIMITCYCLE.ncol,MODELINFO.MATCONTLIMITCYCLE.ncol-1,1)],1);
            trj(counter).x=[ans(:).' 1];
            trj(counter).x0=0;
            trj(counter).timehorizon=xout(MODELINFO.MATCONTLIMITCYCLE.PeriodIdx,ii);
            trj(counter).arcarg=OCMATFINITCONT.arcarg;
            trj(counter).arcinterval=[0 xout(MODELINFO.MATCONTLIMITCYCLE.PeriodIdx,ii)];
            trj(counter).arcposition=[1;ntstcolp1];
            trj(counter).modelparameter=MODELINFO.MATCONTLIMITCYCLE.P0;
            trj(counter).modelparameter(MODELINFO.MATCONTLIMITCYCLE.ActiveParams)=xout(MODELINFO.MATCONTLIMITCYCLE.pars(end),ii);
            trj(counter).modelname=modelnme;
            trj(counter).solverinfo.tmesh=trj(counter).x;
            trj(counter).solverinfo.coeff=xout(:,ii);
            trj(counter).solverinfo.tangent=vout(:,ii);
            if ~isempty(linearization)
                trj(counter).linearization=reshape(linearization(:,ii),MODELINFO.MATCONTLIMITCYCLE.nphase,MODELINFO.MATCONTLIMITCYCLE.nphase);
            else
                trj(counter).linearization=[];
            end
        end
        ocResultStruct.ContinuationSolution=trj;
        if ~counter
            ocResultStruct=[];
            ocResultFieldName='';
            return
        end
        LastodeObj=changeparametervalue(odeObj,trj(counter).modelparameter);
        perStruct.octrajectory=trj(counter);
        perStruct.period=trj(counter).arcinterval(2);

        ocResultStruct.LastLimitCycle=dynprimitive(perStruct,LastodeObj);
        ocResultStruct.LastModel=LastodeObj;
        ocResultStruct.ContinuationInformation=sout;

        ocResultStruct.ContinuationClassification=fieldname;
        ocResultFieldName='MatContContinuation';
    otherwise
        ocResultStruct=fieldvalue;
        ocResultFieldName=fieldname;
end

function [globalvarfile,resultfile,varargout]=loadglobal(odeObj,fieldname)
switch fieldname
    case {'saddlepath2ep'}
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
            resultfile=fullfile(getocmatpath,getocmatfolder(odeObj,'userdata'),['SaveIntermediateResults4' fieldname '.mat']);
        end
        if isempty(OCMATAE) || ~exist(globalvarfile,'file')
            globalvarfile=fullfile(getocmatpath,getocmatfolder(odeObj,'userdata'),['SaveIntermediateResultsGlobalVariable4' fieldname '.mat']);
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
    case {'odesolution','userbvp','odesolution4t','odesolution4s','odesolution4p'}
        global OCMATCONT OCMATODEBVP OCBVP 
        % store value of the global variables OCMATCONT and OCMATODEBVP to
        % temporary variables, in an intermediate step the global variables
        % are changed to values stored in the continuation files and reset
        % afterwards
        if ~isempty(OCMATCONT)
            OCMATCONTORIGINAL=OCMATCONT;
        else
            OCMATCONTORIGINAL=[];
        end
        if ~isempty(OCMATODEBVP)
            OCMATODEBVPORIGINAL=OCMATODEBVP;
        else
            OCMATODEBVPORIGINAL=[];
        end
        if ~isempty(OCBVP)
            OCBVPORIGINAL=OCBVP;
        else
            OCBVPORIGINAL=[];
        end

        % the file names containing continuation information are build by
        % the model specific basic names and the continuation
        % classification, e.g. 4extremal2ep
        if ~isempty(OCMATODEBVP)
            resultfile=[OCMATODEBVP.basicresultfilename fieldname '.mat'];
            globalvarfile=[OCMATODEBVP.basicglobalvarfilename fieldname '.mat'];
        else
            resultfile=fullfile(getocmatpath,getocmatfolder(odeObj,'userdata'),['SaveIntermediateResults4' fieldname '.mat']);
            globalvarfile=fullfile(getocmatpath,getocmatfolder(odeObj,'userdata'),['SaveIntermediateResultsGlobalVariable4' fieldname '.mat']);
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
        % continuation file as MODELINFO.OCMATCONT and MODELINFO.OCMATODEBVP.
        % If the storing process is interrupted before the global variables
        % are reset to its original value the user is informed that s/he
        % has to newly initialize an actual continuation process.
        ocmatmsg('Global variable ''OCMATCONT'' changed.\n')
        OCMATCONT=MODELINFO.OCMATCONT;
        OCMATCONT.meshadaptationflag=1;
        ocmatmsg('Global variable ''OCMATODEBVP'' changed.\n')
        OCMATODEBVP=MODELINFO.OCMATODEBVP;
        ocmatmsg('Global variable ''OCBVP'' changed.\n')
        OCBVP=MODELINFO.OCBVP;
        varargout{1}=OCMATCONTORIGINAL;
        varargout{2}=OCMATODEBVPORIGINAL;
        varargout{3}=OCBVPORIGINAL;
    case {'limitcycle','limitcycleuser','limitcycleuserI'}
        global OCMATCONT OCMATLC OCBVP 
        % store value of the global variables OCMATCONT and OCMATLC to
        % temporary variables, in an intermediate step the global variables
        % are changed to values stored in the continuation files and reset
        % afterwards
        if ~isempty(OCMATCONT)
            OCMATCONTORIGINAL=OCMATCONT;
        else
            OCMATCONTORIGINAL=[];
        end
        if ~isempty(OCMATLC)
            OCMATLCORIGINAL=OCMATLC;
        else
            OCMATLCORIGINAL=[];
        end
        if ~isempty(OCBVP)
            OCBVPORIGINAL=OCBVP;
        else
            OCBVPORIGINAL=[];
        end

        % the file names containing continuation information are build by
        % the model specific basic names and the continuation
        % classification, e.g. 4extremal2ep
        if ~isempty(OCMATLC)
            resultfile=[OCMATLC.basicresultfilename '4' fieldname '.mat'];
            globalvarfile=[OCMATLC.basicglobalvarfilename '4' fieldname '.mat'];
        else
            resultfile=fullfile(getocmatpath,getocmatfolder(odeObj,'userdata'),['SaveIntermediateResults4' fieldname '.mat']);
            globalvarfile=fullfile(getocmatpath,getocmatfolder(odeObj,'userdata'),['SaveIntermediateResultsGlobalVariable4' fieldname '.mat']);
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
        % continuation file as MODELINFO.OCMATCONT and MODELINFO.OCMATLC.
        % If the storing process is interrupted before the global variables
        % are reset to its original value the user is informed that s/he
        % has to newly initialize an actual continuation process.
        ocmatmsg('Global variable ''OCMATCONT'' changed.\n')
        OCMATCONT=MODELINFO.OCMATCONT;
        OCMATCONT.meshadaptationflag=1;
        ocmatmsg('Global variable ''OCMATLC'' changed.\n')
        OCMATLC=MODELINFO.OCMATLC;
        ocmatmsg('Global variable ''OCBVP'' changed.\n')
        OCBVP=MODELINFO.OCBVP;
        varargout{1}=OCMATCONTORIGINAL;
        varargout{2}=OCMATLCORIGINAL;
        varargout{3}=OCBVPORIGINAL;
end

function resetglobal(fieldname,varargin)
switch fieldname
    case {'saddlepath2ep'}
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
    case {'odesolution','userbvp','odesolution4t','odesolution4s','odesolution4p'}
        OCMATCONTORIGINAL=varargin{1};
        OCMATODEBVPORIGINAL=varargin{2};
        OCBVPORIGINAL=varargin{3};
        % reset global variables to its original value or delete it if no
        % global variables existed
        if ~isempty(OCMATCONTORIGINAL)
            ocmatmsg('Global variable ''OCMATCONT'' reset.\n')
            OCMATCONT=OCMATCONTORIGINAL;
        else
            clear global OCMATCONT
        end
        if ~isempty(OCMATODEBVPORIGINAL)
            ocmatmsg('Global variable ''OCMATODEBVP'' reset.\n')
            OCMATODEBVP=OCMATODEBVPORIGINAL;
        else
            clear global OCMATODEBVP
        end
        if ~isempty(OCBVPORIGINAL)
            ocmatmsg('Global variable ''OCBVP'' reset.\n')
            OCBVP=OCBVPORIGINAL;
        else
            clear global OCBVP
        end
        
        
    case {'limitcycle','limitcycleuser','limitcycleuserI'}
        OCMATCONTORIGINAL=varargin{1};
        OCMATLCORIGINAL=varargin{2};
        OCBVPORIGINAL=varargin{3};
        % reset global variables to its original value or delete it if no
        % global variables existed
        if ~isempty(OCMATCONTORIGINAL)
            ocmatmsg('Global variable ''OCMATCONT'' reset.\n')
            OCMATCONT=OCMATCONTORIGINAL;
        else
            clear global OCMATCONT
        end
        if ~isempty(OCMATLCORIGINAL)
            ocmatmsg('Global variable ''OCMATLC'' reset.\n')
            OCMATLC=OCMATLCORIGINAL;
        else
            clear global OCMATLC
        end
        if ~isempty(OCBVPORIGINAL)
            ocmatmsg('Global variable ''OCBVP'' reset.\n')
            OCBVP=OCBVPORIGINAL;
        else
            clear global OCBVP
        end
end

function [solt sole soln sol]=findspecificsolution(sout)
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
