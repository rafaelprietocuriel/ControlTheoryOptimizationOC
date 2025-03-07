function varargout=store(ocObj,varargin)
%
% STORE results of continuation process
%
% STORE(OCOBJ,PROBLEMTYPE) store results of continuation process
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


if isempty(ocObj) || nargin==1
    return
end

fieldname=[];
fieldvalue=[];
if ~ischar(varargin{1}) && nargin==2
    if isdynprimitive(varargin{1}) ||  isocasymptotic(varargin{1}) || isoctrajectory(varargin{1}) || isoccurve(varargin{1}) ||  isocmultipath(varargin{1})
        fieldname=class(varargin{1});
    elseif iscell(varargin{1})
        for ii=1:numel(varargin{1})
            store(ocObj,varargin{1}{ii},varargin{2:end});
        end

        if ~nargout
            assignin('caller',inputname(1),ocObj);
        else
            varargout{1}=ocObj;
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
[ocStruct,ocResultFieldName]=generateelement(ocObj,fieldname,fieldvalue);
if isempty(ocResultFieldName) || isempty(ocStruct)
    return
end
if isfield(ocObj.Result,ocResultFieldName)
    ocObj.Result.(ocResultFieldName){end+1}=ocStruct;
else
    ocObj.Result(1).(ocResultFieldName){1}=ocStruct;
end

if ~nargout
    assignin('caller',inputname(1),ocObj);
else
    varargout{1}=ocObj;
end

%%
function [ocResultStruct,ocResultFieldName]=generateelement(ocObj,fieldname,fieldvalue)
ocResultStruct=[];
ocResultFieldName='';
switch fieldname
%%%%%
%% solutions calculated by GRADIENT method
%%%%%

    case {'staticoptimization','staticindifferencesolution'}
        global OCGRADCONT OCGRADSOL
        [globalvarfile,resultfile,OCGRADCONTORIGINAL,OCGRADSOLORIGINAL]=loadglobal(ocObj,fieldname);

        % load data of information files into workspace
        load(globalvarfile)
        load(resultfile)

        mainfunch=feval(str2func(fieldname));
        % determine if non admissible solution occurs
        [soltarget,sollast,soln,sol]=findspecificsolution(sout);

        counter=0;
        for ii=1:numel(gradout)
            try
                counter=counter+1;
                ocExStruct(counter)=mainfunch{22}(gradout(ii).coeff,gradout(ii).extremal,gradout(ii).tangent);
                ocExStruct(counter).extremal.stepwidth=gradout(ii).stepwidth;
            end
        end
        %OCMATCONT.meshadaptationflag=meshadaptationflag;
        if ~isempty(soltarget)
            sol=soltarget;
            ocResultStruct.ExtremalSolution=ocgradtrajectory(sol.extremal);
            ocResultStruct.NonadmissibleSolution=ocgradtrajectory([]);
        elseif ~isempty(soln)
            ocResultStruct.NonadmissibleSolution=ocgradtrajectory(soln.extremal);
            ocResultStruct.ExtremalSolution=ocgradtrajectory(sol.extremal);
        else
            sol=ocExStruct(end);
            ocResultStruct.ExtremalSolution=ocgradtrajectory(sol.extremal);
            ocResultStruct.NonadmissibleSolution=ocgradtrajectory([]);
        end
        ocResultStruct.ContinuationSolution=ocExStruct;
        ocResultStruct.ContinuationClassification=fieldname;
        ocResultStruct.ContinuationInformation=sout;
        ocResultStruct.GlobalInformation=OCGRADCONT;
        ocResultStruct.ModelInformation=OCGRADSOL;
        ocResultFieldName='Continuation';

        resetglobal(fieldname,OCGRADCONTORIGINAL,OCGRADSOLORIGINAL);

    case {'semismoothoptimization'}
        global OCSTATCONT OCSTATSOL
        [globalvarfile,resultfile,OCSTATCONTORIGINAL,OCSTATSOLORIGINAL]=loadglobal(ocObj,fieldname);

        % load data of information files into workspace
        load(globalvarfile)
        load(resultfile)

        mainfunch=feval(str2func(fieldname));
        % determine if non admissible solution occurs
        [soltarget,sollast,soln,sol]=findspecificsolution(sout);

        counter=0;
        for ii=1:numel(contout)
            try
                counter=counter+1;
                ocExStruct(counter)=mainfunch{22}(contout(ii).coeff,contout(ii).tangent);
            end
        end
        %OCMATCONT.meshadaptationflag=meshadaptationflag;
        if ~isempty(soltarget)
            sol=soltarget;
            ocResultStruct.ExtremalSolution=staticextremal(sol);
            ocResultStruct.NonadmissibleSolution=staticextremal([]);
        elseif ~isempty(soln)
            ocResultStruct.NonadmissibleSolution=staticextremal(soln);
            ocResultStruct.ExtremalSolution=staticextremal(sol);
        else
            sol=ocExStruct(end);
            ocResultStruct.ExtremalSolution=staticextremal(sol);
            ocResultStruct.NonadmissibleSolution=staticextremal([]);
        end
        ocResultStruct.ContinuationSolution=ocExStruct;
        ocResultStruct.ContinuationClassification=fieldname;
        ocResultStruct.ContinuationInformation=sout;
        ocResultStruct.GlobalInformation=OCSTATCONT;
        ocResultStruct.ModelInformation=OCSTATSOL;
        ocResultFieldName='Continuation';

        resetglobal(fieldname,OCSTATCONTORIGINAL,OCSTATSOLORIGINAL);

        
    otherwise
        ocResultStruct=fieldvalue;
        ocResultFieldName=fieldname;
end

function [globalvarfile,resultfile,varargout]=loadglobal(ocObj,fieldname)
switch fieldname
        
%%%%%%%%%%
%% solutions calculated by GRADIENT method
%%%%%%%%%

    case {'staticoptimization','staticindifferencesolution'}
        global OCGRADCONT OCGRADSOL
        % store value of the global variables OCMATCONT and OCMATAE to
        % temporary variables, in an intermediate step the global variables
        % are changed to values stored in the continuation files and reset
        % afterwards
        if ~isempty(OCGRADCONT)
            OCGRADCONTORIGINAL=OCGRADCONT;
        else
            OCGRADCONTORIGINAL=[];
        end
        if ~isempty(OCGRADSOL)
            OCGRADSOLORIGINAL=OCGRADSOL;
        else
            OCGRADSOLORIGINAL=[];
        end

        % the file names containing continuation information are build by
        % the model specific basic names and the continuation
        % classification, e.g. 4extremal2ep
        if ~isempty(OCGRADSOL)
            resultfile=[OCGRADSOL.basicresultfilename '4' fieldname '.mat'];
            globalvarfile=[OCGRADSOL.basicglobalvarfilename '4' fieldname '.mat'];
        else
            resultfile=fullfile(getocmatpath,getocmatfolder(ocObj,'userdata'),['SaveIntermediateResults4' fieldname '.mat']);
            globalvarfile=fullfile(getocmatpath,getocmatfolder(ocObj,'userdata'),['SaveIntermediateResultsGlobalVariable4' fieldname '.mat']);
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
            ocmatmsg('Data file has already been generated at %s.\n',datestr(info.datenum))
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
        ocmatmsg('Global variable ''OCGRADCONT'' changed.\n')
        OCGRADCONT=MODELINFO.OCGRADCONT;
        ocmatmsg('Global variable ''OCGRADSOL'' changed.\n')
        OCGRADSOL=MODELINFO.OCGRADSOL;
        varargout{1}=OCGRADCONTORIGINAL;
        varargout{2}=OCGRADSOLORIGINAL;

        
    case {'semismoothoptimization'}
        global OCSTATCONT OCSTATSOL
        % store value of the global variables OCMATCONT and OCMATAE to
        % temporary variables, in an intermediate step the global variables
        % are changed to values stored in the continuation files and reset
        % afterwards
        if ~isempty(OCSTATCONT)
            OCSTATCONTORIGINAL=OCSTATCONT;
        else
            OCSTATCONTORIGINAL=[];
        end
        if ~isempty(OCSTATSOL)
            OCSTATSOLORIGINAL=OCSTATSOL;
        else
            OCSTATSOLORIGINAL=[];
        end

        % the file names containing continuation information are build by
        % the model specific basic names and the continuation
        % classification, e.g. 4extremal2ep
        if ~isempty(OCSTATSOL)
            resultfile=[OCSTATSOL.basicresultfilename '4' fieldname '.mat'];
            globalvarfile=[OCSTATSOL.basicglobalvarfilename '4' fieldname '.mat'];
        else
            resultfile=fullfile(getocmatpath,getocmatfolder(ocObj,'userdata'),['SaveIntermediateResults4' fieldname '.mat']);
            globalvarfile=fullfile(getocmatpath,getocmatfolder(ocObj,'userdata'),['SaveIntermediateResultsGlobalVariable4' fieldname '.mat']);
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
            ocmatmsg('Data file has already been generated at %s.\n',datestr(info.datenum))
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
        ocmatmsg('Global variable ''OCSTATCONT'' changed.\n')
        OCSTATCONT=MODELINFO.OCSTATCONT;
        ocmatmsg('Global variable ''OCSTATSOL'' changed.\n')
        OCSTATSOL=MODELINFO.OCSTATSOL;
        varargout{1}=OCSTATCONTORIGINAL;
        varargout{2}=OCSTATSOLORIGINAL;

end

function resetglobal(fieldname,varargin)
switch fieldname
        
%%%%%%%%%%
%% solutions calculated by GRADIENT method
%%%%%%%%%

    case {'staticoptimization','staticindifferencesolution'}
        OCGRADCONTORIGINAL=varargin{1};
        OCGRADSOLORIGINAL=varargin{2};
        % reset global variables to its original value or delete it if no
        % global variables existed
        if ~isempty(OCGRADCONTORIGINAL)
            ocmatmsg('Global variable ''OCGRADCONT'' reset.\n')
            OCGRADCONT=OCGRADCONTORIGINAL;
        else
            clear global OCGRADCONT
        end
        if ~isempty(OCGRADSOLORIGINAL)
            ocmatmsg('Global variable ''OCGRADSOL'' reset.\n')
            OCGRADSOL=OCGRADSOLORIGINAL;
        else
            clear global OCGRADSOL
        end

    case {'semismoothoptimization'}
        OCSTATCONTORIGINAL=varargin{1};
        OCSTATSOLORIGINAL=varargin{2};
        % reset global variables to its original value or delete it if no
        % global variables existed
        if ~isempty(OCSTATCONTORIGINAL)
            ocmatmsg('Global variable ''OCSTATCONT'' reset.\n')
            OCSTATCONT=OCSTATCONTORIGINAL;
        else
            clear global OCSTATCONT
        end
        if ~isempty(OCSTATSOLORIGINAL)
            ocmatmsg('Global variable ''OCSTATSOL'' reset.\n')
            OCSTATSOL=OCSTATSOLORIGINAL;
        else
            clear global OCSTATSOL
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
