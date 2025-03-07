function varargout=store(docObj,varargin)
%
% STORE results of continuation process
%
% STORE(OCOBJ,PROBLEMTYPE) store results of continuation process
% given by PROBLEMTYPE:
%   'dextremal2fp'   ... saddle-path of an equilibrium, continuing along the
%                       initial point
%   'dextremalp2ep'  ... saddle-path of an equilibrium, continuing along a
%                       parameter value
%   'dextremalt2ep'  ... saddle-path of an equilibrium, continuing the
%                       truncation time
%   'indifferencesolution'  ... continuation of an indifference threshold
%                       (Skiba point)
%   'limitdextremal' ... continuation of a limitpoint solution
%
% STORE(OCOBJ,OCELEMENT) stores the object OCELEMENT in the Result field of
% the stdocmodel object OCOBJ. 
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
% If no output argument is specified it is written into the object OCOBJ.
%
% OCOBJN=STORE(...) the result is returned to the new stdocmodel instance
% OCOBJN 


if isempty(docObj) || nargin==1
    return
end

fieldname=[];
fieldvalue=[];
if ~ischar(varargin{1}) && nargin==2
    if ismapprimitive(varargin{1}) ||  isdocasymptotic(varargin{1}) || isdoctrajectory(varargin{1}) || isoccurve(varargin{1})
        fieldname=class(varargin{1});
    elseif iscell(varargin{1})
        for ii=1:numel(varargin{1})
            store(docObj,varargin{1}{ii},varargin{2:end});
        end

        if ~nargout
            assignin('caller',inputname(1),docObj);
        else
            varargout{1}=docObj;
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
[ocStruct,ocResultFieldName]=generateelement(docObj,fieldname,fieldvalue);
if isempty(ocResultFieldName)
    return
end
if isfield(docObj.Result,ocResultFieldName)
    docObj.Result.(ocResultFieldName){end+1}=ocStruct;
else
    docObj.Result(1).(ocResultFieldName){1}=ocStruct;
end

if ~nargout
    assignin('caller',inputname(1),docObj);
else
    varargout{1}=docObj;
end

%%
function [ocResultStruct,ocResultFieldName]=generateelement(docObj,fieldname,fieldvalue)
ocResultStruct=[];
ocResultFieldName='';
switch fieldname
    case 'mapprimitive'
        if isfixpoint(fieldvalue)
            ocResultStruct=fieldvalue;
            ocResultFieldName='Fixpoint';
        end
        
    case 'docasymptotic'
        ocResultStruct=fieldvalue;
        ocResultFieldName='ExtremalSolution';
        
    case {'dextremal2fp'}
        global OCMATCONT OCMATAE OCBVP

        [globalvarfile,resultfile,OCMATCONTORIGINAL,OCMATAEORIGINAL,OCBVPORIGINAL]=loadglobal(docObj,fieldname);

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
                mainfunch{6}(bvpout(ii).tmesh,bvpout(ii).coeff,bvpout(ii).tangent);
                ocExStruct(counter)=mainfunch{22}(bvpout(ii).tmesh,bvpout(ii).coeff,bvpout(ii).tangent);
            end
        end
        ocFP=mapprimitive(OCMATAE.saddlepoint,OCMATCONT.HE.arcarg(end),OCMATAE.linearization);
        if ~isempty(soln)
            mainfunch{6}(soln.solverinfo.tmesh,soln.solverinfo.coeff,soln.solverinfo.tangent);
            soln=mainfunch{22}(soln.solverinfo.tmesh,soln.solverinfo.coeff,soln.solverinfo.tangent);
            mainfunch{6}(sol.solverinfo.tmesh,sol.solverinfo.coeff,sol.solverinfo.tangent);
            sol=mainfunch{22}(soln.solverinfo.tmesh,soln.solverinfo.coeff,soln.solverinfo.tangent);
            ocResultStruct.NonadmissibleSolution=docasymptotic(doctrajectory(soln,docObj),ocFP);
            ocResultStruct.ExtremalSolution=docasymptotic(doctrajectory(sol,docObj),ocFP);
        else
            if ~isempty(soltarget)
                mainfunch{6}(soltarget.solverinfo.tmesh,soltarget.solverinfo.coeff,soltarget.solverinfo.tangent);
                soltarget=mainfunch{22}(soltarget.solverinfo.tmesh,soltarget.solverinfo.coeff,soltarget.solverinfo.tangent);
                ocResultStruct.ExtremalSolution=docasymptotic(doctrajectory(soltarget,docObj),ocFP);
            else
                ocResultStruct.ExtremalSolution=docasymptotic(doctrajectory(ocExStruct(end),docObj),ocFP);
            end
            ocResultStruct.NonadmissibleSolution=docasymptotic([]);
        end
        ocResultStruct.ContinuationSolution=ocExStruct;
        ocResultStruct.ContinuationClassification=fieldname;
        ocResultStruct.ContinuationInformation=sout;
        ocResultStruct.LimitSet=ocFP;
        ocResultFieldName='Continuation';

        resetglobal(fieldname,OCMATCONTORIGINAL,OCMATAEORIGINAL,OCBVPORIGINAL);

    case 'dextremalp2fp'
        global OCMATCONT OCMATAE OCBVP
        [globalvarfile,resultfile,OCMATCONTORIGINAL,OCMATAEORIGINAL,OCBVPORIGINAL]=loadglobal(docObj,fieldname);

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
                mainfunch{6}(bvpout(ii).tmesh,bvpout(ii).coeff,bvpout(ii).tangent);
                ocExStruct(counter)=mainfunch{22}(bvpout(ii).tmesh,bvpout(ii).coeff,bvpout(ii).tangent);
            end
        end

        if ~isempty(soln)
            equcoord=soln.solverinfo.parametercoord(soln.solverinfo.equilibriumcoord);
            docObjtmp=changeparametervalue(docObj,OCMATAE.varyparameterindex,soln.solverinfo.continuationparameter);
            hatx.y=soln.solverinfo.coeff(equcoord);
            hatx.x=0;
            hatx.arcarg=OCMATCONT.HE.arcarg(end);
            docFP=mapprimitive(hatx);
            hatx.linearization=jacobian(docObjtmp,docFP);
            docFP=mapprimitive(hatx);
            ocResultStruct.NonadmissibleSolution=docasymptotic(doctrajectory(soln,docObjtmp),docFP);
            equcoord=sol.solverinfo.parametercoord(sol.solverinfo.equilibriumcoord);
            docObjtmp=changeparametervalue(docObj,OCMATAE.varyparameterindex,sol.solverinfo.continuationparameter);
            hatx.y=sol.solverinfo.coeff(equcoord);
            hatx.x=0;
            hatx.arcarg=OCMATCONT.HE.arcarg(end);
            docFP=mapprimitive(hatx);
            hatx.linearization=jacobian(docObjtmp,docFP);
            docFP=mapprimitive(hatx);
            ocResultStruct.ExtremalSolution=docasymptotic(doctrajectory(sol,docObjtmp),docFP);
        else
            if ~isempty(soltarget)
                sol=soltarget;
            else
                sol=ocExStruct(end);
            end
            equcoord=sol.solverinfo.parametercoord(sol.solverinfo.equilibriumcoord);
            docObjtmp=changeparametervalue(docObj,OCMATAE.varyparameterindex,sol.solverinfo.continuationparameter);
            hatx.y=sol.solverinfo.coeff(equcoord);
            hatx.x=0;
            hatx.arcarg=OCMATCONT.HE.arcarg(end);
            docFP=mapprimitive(hatx);
            hatx.linearization=jacobian(docObjtmp,docFP);
            docFP=mapprimitive(hatx);
            ocResultStruct.ExtremalSolution=docasymptotic(doctrajectory(sol,docObj),docFP);
            ocResultStruct.NonadmissibleSolution=ocasymptotic([]);
        end
        ocResultStruct.ContinuationParameter=OCMATAE.varyparameterindex;
        ocResultStruct.ContinuationSolution=ocExStruct;
        ocResultStruct.ContinuationClassification=fieldname;
        ocResultStruct.ContinuationInformation=sout;
        ocResultStruct.LastModel=docObjtmp;
        ocResultFieldName='Continuation';

        resetglobal(fieldname,OCMATCONTORIGINAL,OCMATAEORIGINAL,OCBVPORIGINAL);
        
    case {'dindifferencesolution'}
        global OCMATCONT OCMATINDIF OCBVP

        [globalvarfile,resultfile,OCMATCONTORIGINAL,OCMATINDIFORIGINAL,OCBVPORIGINAL]=loadglobal(docObj,fieldname);

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
                mainfunch{6}(bvpout(ii).tmesh,bvpout(ii).coeff,bvpout(ii).tangent);
                ocExStruct(counter)=mainfunch{22}(bvpout(ii).tmesh,bvpout(ii).coeff,bvpout(ii).tangent);
                indifferencepointStruct.y(1:OCBVP.nummap,counter)=ocExStruct(ii).y0(1:OCBVP.nummap);
                indifferencepointStruct.arcarg(counter)=OCMATCONT.HE.arcarg(1);
                indifferencepointStruct.userinfo.continuationparameter(counter)=bvpout(ii).coeff(end);
            end
        end
        indifferencepointStruct.arcposition=[];
        for ii=1:OCMATINDIF.indifferenceorder
            ocFP{ii}=mapprimitive(OCMATINDIF.saddlepoint{ii},OCMATINDIF.arcarg(OCMATINDIF.arcargcoord(2,ii)),OCMATINDIF.linearization{ii});
        end

        if ~isempty(soln)
            ocResultStruct.NonadmissibleSolution=sol2docmultipath(soln,ocFP);
            ocResultStruct.ExtremalSolution=sol2docmultipath(sol,ocFP);
        else
            if ~isempty(soltarget)
                ocResultStruct.ExtremalSolution=sol2docmultipath(soltarget,ocFP);
            else
                ocResultStruct.ExtremalSolution=sol2docmultipath(ocExStruct(end),ocFP);
            end
            ocResultStruct.NonadmissibleSolution=docasymptotic([]);
        end
        if strcmp(fieldname,'dindifferencesolution')
            ocResultStruct.IndifferencePointCurve=occurve(indifferencepointStruct);
        end
        ocResultStruct.ContinuationSolution=ocExStruct;
        ocResultStruct.ContinuationClassification=fieldname;
        ocResultStruct.ContinuationInformation=sout;
        ocResultFieldName='Continuation';

        resetglobal(fieldname,OCMATCONTORIGINAL,OCMATINDIFORIGINAL,OCBVPORIGINAL);

    otherwise
        ocResultStruct=fieldvalue;
        ocResultFieldName=fieldname;
end

function [globalvarfile,resultfile,varargout]=loadglobal(docObj,fieldname)
switch fieldname
    case {'dextremal2fp','dextremalp2fp'}
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
        % classification, e.g. 4dextremal2ep
        if ~isempty(OCMATAE)
            resultfile=[OCMATAE.basicresultfilename '4' fieldname '.mat'];
            globalvarfile=[OCMATAE.basicglobalvarfilename '4' fieldname '.mat'];
        end

        if isempty(OCMATAE) || ~exist(resultfile,'file')
            resultfile=fullfile(getocmatpath,getocmatfolder(docObj,'userdata'),['SaveIntermediateResults4' fieldname '.mat']);
        end
        if isempty(OCMATAE) || ~exist(globalvarfile,'file')
            globalvarfile=fullfile(getocmatpath,getocmatfolder(docObj,'userdata'),['SaveIntermediateResultsGlobalVariable4' fieldname '.mat']);
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
        
    case {'dindifferencesolution'}
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
            OCMATINDIF=[];
        end

        % the file names containing continuation information are build by
        % the model specific basic names and the continuation
        % classification, e.g. 4dextremal2ep
        if ~isempty(OCMATINDIF)
            resultfile=[OCMATINDIF.basicresultfilename '4' fieldname '.mat'];
            globalvarfile=[OCMATINDIF.basicglobalvarfilename '4' fieldname '.mat'];
        else
            resultfile=fullfile(getocmatpath,getocmatfolder(docObj,'userdata'),['SaveIntermediateResults4' fieldname '.mat']);
            globalvarfile=fullfile(getocmatpath,getocmatfolder(docObj,'userdata'),['SaveIntermediateResultsGlobalVariable4' fieldname '.mat']);
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
    case {'dextremal2fp','dextremalp2fp'}
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
        
    case {'dindifferencesolution'}
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
