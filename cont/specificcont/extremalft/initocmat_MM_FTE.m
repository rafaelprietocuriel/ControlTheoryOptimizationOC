function sol=initocmat_MM_FTE(mmObj,ocTrj,continuationtype,continuationindex,continuationtarget,varargin)
%
% initocmat_FTE_T initialization for asymptotic extremal calculation
%
% SOL=initocmat_FTE_T(OCOBJ,OCASYM,CONTCOORDINATE,TARGETVALUE) the
% continuation process is initialized (started) at a solution OCASYM
% (instance of the ocasymptotic class) that (usually) has been calculated
% during a previous continuation. Therefore, neither the 'TruncationTime' ,
% nor the 'PathType' (see initocmat_FTE_T) can be provided. This
% information is taken from the OCASYM object. 
%
% For a more detailed discussion see INITOCMAT_AE_EP
%
clear global OCMATCONT OCMATFTE
global OCMATCONT OCMATFTE
sol=[];
initialstate=[];
initialcostate=[];
initialarcargument=[];
findoptimalswitchingtime=[];
fixendstate=[];
fixinitialstate=[];
freeparametervector=[];
findoptimalparameter=[];
exogenousfunction=[];
connectiontimeindex=[];
userbc=[];
hitvalue=[];

objectivevaluecalc=[];
explicitconnectime=[];
addpart=[];
opt=[];
if isempty(mmObj)
    ocmatmsg('oc model is empty.')
    return
end

% contiuationtype: time, state, parameter
% continuationindex: (time) number of switching/endtime, (state)
% coordinate, (parameter) index or name of parameter
% continuationtarget: value
optionidx=find(strcmpi(varargin,'option'));
addpartidx=find(strcmpi(varargin,'addpart'));
initialstateidx=find(strcmpi(varargin,'initialstate'));
initialcostateidx=find(strcmpi(varargin,'initialcostate'));
initialarcargumentidx=find(strcmpi(varargin,'initialarcargument'));
fixinitialstateidx=find(strcmpi(varargin,'fixinitialstate'));
fixendstateidx=find(strcmpi(varargin,'fixendstate'));
freeparametervectoridx=find(strcmpi(varargin,'freeparametervector'));
findoptimalswitchingtimeidx=find(strcmpi(varargin,'findoptimalswitchingtime'));
findoptimalparameteridx=find(strcmpi(varargin,'findoptimalparameter'));
objectivevaluecalcidx=find(strcmpi(varargin,'objectivevaluecalc'));
connectiontimeindexidx=find(strcmpi(varargin,'connectiontimeindex'));
explicitconnectimeidx=find(strcmpi(varargin,'explicitconnectime'));
userbcidx=find(strcmpi(varargin,'userbc'));
hitvalueidx=find(strcmpi(varargin,'hitvalue'));
exogenousfunctionidx=find(strcmpi(varargin,'exogenousfunction'));
exogenousinitialstatesidx=find(strcmpi(varargin,'exogenousinitialstates'));
if ~isempty(addpartidx)
    addpart=varargin{addpartidx+1};
end
if ~isempty(fixendstateidx)
    fixendstate=varargin{fixendstateidx+1};
end
if ~isempty(fixinitialstateidx)
    fixinitialstate=varargin{fixinitialstateidx+1};
end
if ~isempty(freeparametervectoridx)
    freeparametervector=varargin{freeparametervectoridx+1};
end
if ~isempty(findoptimalparameteridx)
    findoptimalparameter=varargin{findoptimalparameteridx+1};
end
if ~isempty(objectivevaluecalcidx)
    objectivevaluecalc=varargin{objectivevaluecalcidx+1};
end
if ~isempty(explicitconnectimeidx)
    explicitconnectime=varargin{explicitconnectimeidx+1};
end
if ~isempty(connectiontimeindexidx)
    connectiontimeindex=varargin{connectiontimeindexidx+1};
end
if ~isempty(findoptimalswitchingtimeidx)
    findoptimalswitchingtime=varargin{findoptimalswitchingtimeidx+1};
end
if ~isempty(exogenousfunctionidx)
    exogenousfunction=varargin{exogenousfunctionidx+1};
end
if ~isempty(exogenousinitialstatesidx)
    exogenousinitialstates=varargin{exogenousinitialstatesidx+1};
end
if ~isempty(userbcidx)
    userbc=varargin{userbcidx+1};
end
if ~isempty(hitvalueidx)
    hitvalue=varargin{hitvalueidx+1};
end
if isempty(objectivevaluecalc)
    objectivevaluecalc=0;
end
if isempty(findoptimalswitchingtime)
    findoptimalswitchingtime=false;
end
if isempty(findoptimalparameter)
    findoptimalparameter=false;
end
if addpart
    if isempty(ocTrj)
        numberofstages=1;
    else
        numberofstages=numberofparts(ocTrj)+1;
    end
else
    numberofstages=numberofparts(ocTrj);
end

if numberofstages>numberofmodels(mmObj);
    sol=[];
    return
end
OCMATFTE.numberofstages=numberofstages;
if ~iscell(ocTrj)
    OCMATFTE.optimalswitchingtime=zeros(1,numberofstages);
    if isoctrajectory(ocTrj)
        ocTrjCell{1}=ocTrj;
    elseif ismmultipath(ocTrj)
        for ii=1:numberofparts(ocTrj)
            ocTrjCell{ii}=ocTrj(ii);
        end
        if numberofstages==numberofparts(ocTrj)
            OCMATFTE.optimalswitchingtime=optimaltransition(ocTrj);
        else
            OCMATFTE.optimalswitchingtime=zeros(1,numberofstages);
        end
    end
else
    ocTrjCell=ocTrj;
    OCMATFTE.optimalswitchingtime=zeros(1,numberofstages);
end
switch continuationtype
    case {'initialstate','endstate'}
        OCMATFTE.continuationtype=0;
        
    case {'switchingtime','time'}
        OCMATFTE.continuationtype=1;
    case 'parameter'
        OCMATFTE.continuationtype=2;
        parametervaluetmp=parametervalue(mmObj);
        %         parindex=continuationindex;
        %         if ischar(parindex)
        %             parindex=repmat({parindex},1,numberofstages);
        %         end
        %
        %         if iscell(parindex)
        %             % character: same parameter variables for each stage
        %             % cell:
        %             if iscell(parindex)
        %                 if length(parindex)~=numberofstages
        %                     return
        %                 else
        %                 end
        %             end
        %             continuationindex=parameterindex(mmObj,parindex);
        %             if length(continuationindex{1})==size(continuationtarget,1) && size(continuationtarget,2)==1
        %                 continuationtarget=num2cell(repmat(continuationtarget,1,numberofstages));
        %             elseif size(continuationtarget,2)~=numberofstages
        %                 return
        %             end
        %
        %         else
        %             ocmatmsg('Input argument for parameter has to be a character or cell array.')
        %             return
        %         end
        continuationindextmp=[];
        continuationtargettmp=[];
        if iscell(continuationindex)
            if numberofstages==length(continuationindex)
                continuationindextmp=continuationindex;
                for jj=1:length(continuationindex)
                    if isempty(continuationindex{jj})
                        continuationtargettmp{jj}=[];
                    else
                        if iscell(continuationtarget)
                            continuationtargettmp{jj}=continuationtarget{jj};
                        else
                            if length(continuationtarget)==length(continuationindex)
                                continuationtargettmp{jj}=continuationtarget(jj);
                            else
                                continuationtargettmp{jj}=continuationtarget;
                            end
                        end
                    end
                end
            else
                wrongcellsize=1;
            end
        elseif ischar(continuationindex)
            if ~length(continuationindex)==length(continuationtarget)
                wrongcontinuationtargetstruct=1;
            end
            % it is assumed that for each model the same parameters are free
            for jj=1:numberofstages
                continuationindextmp{jj}=continuationindex;
                continuationtargettmp{jj}=continuationtarget;
            end
        else
            ocmatmsg('Parameter vector has to be a cell-array or character.')
            return
        end
        continuationindex=continuationindextmp;
        continuationtarget=continuationtargettmp;
        counter=0;
        for jj=1:numberofstages
            counter=counter+1;
            if ischar(continuationindex{jj})
                continuationindex{jj}=parameterindex(mmObj(jj),continuationindex{jj});
                OCMATFTE.continuationindex{jj}=continuationindex{jj};
            else
                OCMATFTE.continuationindex{jj}=continuationindex{jj};
            end
            if length(continuationtarget{jj})==length(continuationindex{jj})
                OCMATFTE.initialparameter{jj}=parametervaluetmp{jj}(continuationindex{jj});
                if isempty(continuationtarget{jj})
                    OCMATFTE.continuationvector{jj}=[];
                else
                    OCMATFTE.continuationvector{jj}=continuationtarget{jj}-OCMATFTE.initialparameter{jj};
                end
            else
                ocmatmsg('Target parameter vector cell has wrong size.')
                return
            end
        end
        continuationparameter=0;
end
if ~isempty(freeparametervector)
    nummodel=numberofmodels(mmObj);
    for jj=1:nummodel
        if ischar(freeparametervector{jj})
            OCMATFTE.freeparameterindex{jj}=parameterindex(mmObj(jj),freeparametervector{jj});
            freeparametername{jj}=parametername(mmObj(jj),OCMATFTE.freeparameterindex{jj});
            if ~isempty(freeparametervector{jj})
                freeparametervalue{jj}=parametervalue(mmObj(jj),OCMATFTE.freeparameterindex{jj});
            else
                freeparametervalue{jj}=[];
            end
        else
            OCMATFTE.freeparameterindex{counter}=freeparametervector{jj};
            freeparametervalue{jj}=parametervalue(mmObj(jj),freeparametervector{jj});
        end
    end
    numberoffreeparameter=-inf;
    freeparameternametotal='';
    freeparametervaluetotal=[];
    for jj=1:nummodel
        numberoffreeparameter=max(numberoffreeparameter,length(OCMATFTE.freeparameterindex{jj}));
        freeparameternametotal=[freeparameternametotal;freeparametername{jj}];
        freeparametervaluetotal=[freeparametervaluetotal;freeparametervalue{jj}.'];
    end
    [freeparameternametotal,idx]=unique(freeparameternametotal);
    freeparametervalue=freeparametervaluetotal(idx);
    freeparametervalue=freeparametervalue(:).';
    for jj=1:nummodel
        if ~isempty(freeparametervector{jj})
            actualnames=parametername(mmObj(jj),OCMATFTE.freeparameterindex{jj});
            counter2=0;
            for kk=1:length(actualnames)
                idx=strmatch(actualnames{kk},freeparameternametotal);
                if ~isempty(idx)
                    counter2=counter2+1;
                    freeparametervectorcoordinate{jj}(counter2)=idx;
                end
            end
        else
            freeparametervectorcoordinate{jj}=[];
        end
    end

end
OCMATFTE.freeparametervector=freeparametervector;
OCMATFTE.findoptimalparameter=findoptimalparameter;
OCMATFTE.userbc=userbc;
if ~isempty(hitvalue)
    %OCMATFTE.targetdistance=norm(OCMATFTE.startvalue-hitvalue)/norm(OCMATFTE.startvalue-targetvalue);
    OCMATFTE.hitvalue=hitvalue;
else
     OCMATFTE.targetdistance=[];
     OCMATFTE.hitvalue=[];
end

%%% generate initial octrajectory if not provided
if addpart
    if isempty(ocTrj)|| ((length(ocTrjCell)==numberofstages-1 || (length(ocTrjCell)==numberofstages && isempty(ocTrjCell{numberofstages}))))
        if ~isempty(initialstateidx)
            initialstate=varargin{initialstateidx+1};
            initialstate=initialstate(:);
        end
        if ~isempty(initialcostateidx)
            initialcostate=varargin{initialcostateidx+1};
            initialcostate=initialcostate(:);
        end
        if ~isempty(initialarcargumentidx)
            initialarcargument=varargin{initialarcargumentidx+1};
        end
        if ~isempty(optionidx)
            opt=varargin{optionidx+1};
        end
        if isempty(opt)
            opt=defaultocoptions;
        end
        if isempty(initialarcargument)
            initialarcargument=0;
        end
        if isempty(initialstate)
            initialstate=state(mmObj(numberofstages-1),ocTrjCell{numberofstages-1},1);
            initialstate=initialstate(:,end);
        end
        if isempty(initialcostate) && ~isempty(ocTrj)
            initialcostate=costate(mmObj(numberofstages-1),ocTrjCell{numberofstages-1},1);
            initialcostate=initialcostate(:,end);
        elseif isempty(initialcostate) && isempty(ocTrj)
            initocPt.y=initialstate;
            initocPt.x=0;
            initocPt.arcarg=initialarcargument;
            initocPt.arcinterval=[0 0];
            initocPt.arcposition=[1;1];
            initialcostate=transversalitycondition(mmObj(numberofstages),octrajectory(initocPt));
        end
        n=getocoptions(opt,'GENERAL','TrivialArcMeshNum');
        ocTrjN.x=linspace(0,1,n);
        if numberofstages>1
            t=time(mmObj(numberofstages-1),ocTrjCell{numberofstages-1},1);
        else
            t=0;
        end
        ocTrjN.y=[initialstate(:);initialcostate];
        ocTrjN.y=ocTrjN.y(:,ones(1,n));
        ocTrjN.arcposition=[1;n];
        ocTrjN.arcinterval=[t(end) t(end)];
        ocTrjN.arcarg=initialarcargument;
        ocTrjN.x0=t(end);
        ocTrjN.timehorizon=t(end);
        ocTrjN.modelparameter=parametervalue(mmObj(numberofstages));
        ocTrjN.modelname=modelname(mmObj(numberofstages));
        ocTrjCell{numberofstages}=octrajectory(ocTrjN);
    end
end

if ~isempty(ocTrjCell)
    ocTrjMP=mmultipath(ocTrjCell);
else
    ocTrjMP=ocTrj;
end
OCMATFTE.objectivevaluecalc=objectivevaluecalc;

OCMATFTE.exogenousfunction=exogenousfunction;
if exogenousfunction
    OCMATFTE.exogenousinitialstates=exogenousinitialstates;
    OCMATFTE.exogenousnumberofstates=length(exogenousinitialstates);
end

scoord=statecoord(mmObj);
cscoord=costatecoord(mmObj);
for ii=1:numberofstages
    OCMATFTE.statecostatecoord{ii}=[scoord{ii}(:).' cscoord{ii}(:).'];
end

% dermine if the model explicitly depends on the connection times, in that
% case the

OCMATFTE.dhamiltoniandctcoord=[];
if isempty(explicitconnectime) || explicitconnectime
    OCMATFTE.explicitconnectiontime=any(explicitconnectiontime(mmObj));
else
    OCMATFTE.explicitconnectiontime=explicitconnectiontime;
    OCMATFTE.dhamiltoniandctcoord=[];
end

numberofodes=cscoord{1}(end);
if OCMATFTE.explicitconnectiontime
    OCMATFTE.connectiontimenumber=connectiontimenumber(mmObj);
    OCMATFTE.dhamiltoniandctcoord=2*statenum(mmObj(1))+(1:max(OCMATFTE.connectiontimenumber)); %it is assumed that every part has the same number of states and costates
    OCMATFTE.connectiontimenumber=max(OCMATFTE.connectiontimenumber);
    dHdctval=zeros(OCMATFTE.connectiontimenumber,1);
    dSdctval=zeros(OCMATFTE.connectiontimenumber,1);
    for ii=1:numberofstages
        dSdctval=dSdctval+dsalvagedct(mmObj(ii),ocTrjMP);
    end
    for ii=1:numberofstages
        ocTrjStruct=struct(ocTrjMP(ii));
        if length(ocTrjMP(ii).y(:,1))==2*statenum(mmObj(ii))
            val=dhamiltoniandct(mmObj(ii),ocTrjMP,1);
            if ii==1
                ocTrjStruct.y(OCMATFTE.dhamiltoniandctcoord,:)=repmat(dSdctval,1,size(val,2))+[zeros(OCMATFTE.connectiontimenumber,1) cumsum(val(:,1:end-1)+val(:,2:end),1)/2.*repmat(diff(time(mmObj(ii),ocTrjMP,1)),OCMATFTE.connectiontimenumber,1)];
            else
                ocTrjStruct.y(OCMATFTE.dhamiltoniandctcoord,:)=repmat(dHdctval,1,size(val,2))+[zeros(OCMATFTE.connectiontimenumber,1) cumsum(val(:,1:end-1)+val(:,2:end),1)/2.*repmat(diff(time(mmObj(ii),ocTrjMP,1)),OCMATFTE.connectiontimenumber,1)];
            end
            ocTrjCell{ii}=octrajectory(ocTrjStruct);
            dHdctval=ocTrjStruct.y(OCMATFTE.dhamiltoniandctcoord,end);
        else
            ocTrjCell{ii}=ocTrjMP(ii);
            dHdctval=ocTrjMP(ii).y(OCMATFTE.dhamiltoniandctcoord,end);
        end
    end
    numberofodes=numberofodes+OCMATFTE.connectiontimenumber;
    ocTrjMP=mmultipath(ocTrjCell);
end

if objectivevaluecalc
    Oval=0;
    OCMATFTE.objectivevaluecoord=2*statenum(mmObj(1))+length(OCMATFTE.dhamiltoniandctcoord)+1; %it is assumed that every part has the same number of states and costates
    for ii=1:numberofstages
        ocTrjStruct=struct(ocTrjMP(ii));
        if length(ocTrjMP(ii).y(:,1))<OCMATFTE.objectivevaluecoord
            o=objectivefunction(mmObj(ii),ocTrjMP,1);
            ocTrjStruct.y(OCMATFTE.objectivevaluecoord,:)=Oval+[0 cumsum((o(1:end-1)+o(2:end))/2.*diff(time(mmObj(ii),ocTrjMP,1)))];
            ocTrjCell{ii}=octrajectory(ocTrjStruct);
            Oval=ocTrjStruct.y(OCMATFTE.objectivevaluecoord,end);
        else
            ocTrjCell{ii}=ocTrjMP(ii);
            Oval=ocTrjMP(ii).y(OCMATFTE.objectivevaluecoord,end);
        end
    end
    ocTrjMP=mmultipath(ocTrjCell);
    numberofodes=numberofodes+1;
    
else
    OCMATFTE.objectivevaluecoord=[];
end
if OCMATFTE.exogenousfunction
    OCMATFTE.exogenousdynamicscoordinate=numberofodes+1:numberofodes+OCMATFTE.exogenousnumberofstates;
else
    OCMATFTE.exogenousdynamicscoordinate=[];
end

% initialize global variable (OCMATCONT) for general continuation process
OCMATCONT.modelname=modelname(mmObj);

% initialize global variable (OCMATFTE) for specific continuation process 
% initialize functions 
% empty function handles mean either that they have to be computed
% numerically (derivatives, jacobian, hessian) or that standard functions are
% used (e.g. asymptotic boundary condition, initial condition)
% mode and path specific variables
OCMATFTE.parametervalue=parametervalue(mmObj);
OCMATFTE.autonomous=isautonomous(mmObj);

for ii=1:numberofstages
    OCMATCONT.modelfunc{ii}=modelspecificfunc(mmObj(ii),'4FiniteHorizonPathContinuation');
    funch{ii}=OCMATCONT.modelfunc{ii}(); % model specific function handles for saddle path continuation
end
for ii=1:numberofstages
    if OCMATFTE.explicitconnectiontime
        OCMATFTE.canonicalsystem{ii}=funch{ii}{1}{1};
        OCMATFTE.dhamiltoniandct{ii}=funch{ii}{1}{2};
        OCMATFTE.dsalvagedct{ii}=funch{ii}{1}{3};
    else
        OCMATFTE.canonicalsystem{ii}=funch{ii}{1}{1};
    end
    if OCMATFTE.exogenousfunction
        OCMATFTE.exogenousdynamics{ii}=funch{ii}{1}{4};
        OCMATFTE.exogenousjacobian{ii}=funch{ii}{2}{8};
        OCMATFTE.exogenousparameterjacobian{ii}=funch{ii}{2}{9};
        OCMATFTE.exogenousdct{ii}=funch{ii}{2}{10};
    end
    OCMATFTE.canonicalsystemjacobian{ii}=funch{ii}{2}{1};
    OCMATFTE.canonicalsystemparameterjacobian{ii}=funch{ii}{2}{2};
    if OCMATFTE.explicitconnectiontime
        OCMATFTE.derivativeconnectiontime{ii}=funch{ii}{2}{4};
        OCMATFTE.dhamiltoniandctjacobian{ii}=funch{ii}{2}{5};
        OCMATFTE.dhamiltoniandctparameterjacobian{ii}=funch{ii}{2}{6};
        OCMATFTE.d2hamiltoniandct2{ii}=funch{ii}{2}{7};
    end
    if ~isempty(hitvalue)
        try
            OCMATFTE.hitvaluefunc=funch{ii}{5}{4};
        catch
            OCMATFTE.hitvaluefunc{ii}=[];
        end
    end
    if OCMATFTE.userbc
        OCMATFTE.userfunctionbc{ii}=funch{ii}{5}{9};
    end
    OCMATFTE.canonicalsystemhessian{ii}=funch{ii}{3}{1};
    OCMATFTE.canonicalsystemparameterhessian{ii}=funch{ii}{3}{2};

    % function for the boundary conditions
    OCMATFTE.bcinitial{ii}=funch{ii}{5}{1};
    OCMATFTE.bctransversality{ii}=funch{ii}{5}{2};
    OCMATFTE.bcoptimalhorizon{ii}=funch{ii}{5}{3};
    OCMATFTE.salvagevalue{ii}=funch{ii}{5}{6};
    OCMATFTE.bcconnectingparts{ii}=funch{ii}{5}{7};
    OCMATFTE.bcoptimalconnectingparts{ii}=funch{ii}{5}{8};
    
    % function for Jacobian
    OCMATFTE.bcjacobianinitial{ii}=funch{ii}{6}{1};
    OCMATFTE.bcjacobiantransversality{ii}=funch{ii}{6}{2};

    % function describing the hybrid structure of the problem
    OCMATFTE.hybridinfo{ii}=funch{ii}{7}{1};
    OCMATFTE.domain{ii}=funch{ii}{7}{2};
    OCMATFTE.guard{ii}=funch{ii}{7}{3};
    OCMATFTE.reset{ii}=funch{ii}{7}{4};
    OCMATFTE.switchtime{ii}=funch{ii}{7}{5};
    OCMATFTE.jacobianguard{ii}=funch{ii}{7}{7};
    OCMATFTE.jacobianreset{ii}=funch{ii}{7}{8};
    OCMATFTE.domaindiscretization{ii}=funch{ii}{7}{9};
    
    OCMATFTE.objectivefunction{ii}=funch{ii}{8}{1};
    OCMATFTE.objectivefunctionjacobian{ii}=funch{ii}{8}{2};
    OCMATFTE.objectivefunctionparameterjacobian{ii}=funch{ii}{8}{3};
    OCMATFTE.objectivefunctionderivativetime{ii}=funch{ii}{8}{4};
    if OCMATFTE.explicitconnectiontime
        OCMATFTE.objectivefunctionderivativeconnectiontime{ii}=funch{ii}{8}{5};
    end
    if ~isautonomous(mmObj(ii))
        OCMATFTE.canonicalsystemderivativetime{ii}=funch{ii}{2}{3};
    end
    % general function
    OCMATFTE.plotcontinuation{ii}=funch{ii}{11};
    OCMATFTE.testadmissibility{ii}=funch{ii}{12};
    OCMATFTE.datapath{ii}=funch{ii}{20};
    OCMATFTE.saveintermediatefiles{ii}=funch{ii}{21};

    hybridinfo{ii}=OCMATFTE.hybridinfo{ii}();
    
    for jj=1:numel(hybridinfo{ii}.arcarg)
        domaindata{ii}(jj)=OCMATFTE.domain{ii}(hybridinfo{ii}.arcarg(jj));
    end
    for jj=1:numel(domaindata{ii})
        OCMATCONT.DOMAINDDATA{ii}(jj).numode=domaindata{ii}(jj).odedim;
        OCMATCONT.DOMAINDDATA{ii}(jj).numae=domaindata{ii}(jj).aedim;
        OCMATCONT.DOMAINDDATA{ii}(jj).daeorder=domaindata{ii}(jj).daeorder;
        OCMATCONT.DOMAINDDATA{ii}(jj).numeq=numel(domaindata{ii}(jj).daeorder);%number of equations
        if objectivevaluecalc
            OCMATCONT.DOMAINDDATA{ii}(jj).numode=domaindata{ii}(1).odedim+1;
            if jj==1
                OCMATCONT.DOMAINDDATA{ii}(jj).numeq=OCMATCONT.DOMAINDDATA{ii}(1).numeq+1;%number of equations
            else
                OCMATCONT.DOMAINDDATA{ii}(jj).numeq=OCMATCONT.DOMAINDDATA{ii}(1).numeq;%number of equations
            end
        end
        OCMATCONT.DOMAINDDATA{ii}(jj).eqcoord=1:OCMATCONT.DOMAINDDATA{ii}(jj).numeq;
        OCMATCONT.DOMAINDDATA{ii}(jj).odecoord=1:domaindata{ii}(jj).odedim;
        OCMATCONT.DOMAINDDATA{ii}(jj).aecoord=domaindata{ii}(jj).odedim+(1:domaindata{ii}(jj).aedim);
    end
end

OCMATFTE.findoptimalswitchingtime=findoptimalswitchingtime;

sol=generatesolstruct(ocTrjMP);
if ~objectivevaluecalc && ~OCMATFTE.explicitconnectiontime
    sol.y=sol.y(OCMATFTE.statecostatecoord{1},:);
end

if ~isempty(freeparametervector)
    OCMATFTE.numberoffreeparameter=numberoffreeparameter;
    %OCMATFTE.freeparametervectorcoordinate=length(sol.parameters)+(1:numberoffreeparameter);
    for ii=1:nummodel
       OCMATFTE.freeparametervectorcoordinate{ii}=length(sol.parameters)+freeparametervectorcoordinate{ii};
    end
    
    sol.parameters=[sol.parameters freeparametervalue];
end

depvarFirst=dependentvar(ocTrjMP(1));
depvarLast=dependentvar(ocTrjMP(numberofstages));
OCMATFTE.initialstate=depvarFirst(scoord{1},1);
OCMATFTE.endstate=depvarLast(scoord{numberofstages},end);


arcintv=arcinterval(ocTrjMP);
timepoints=repmat(NaN,1,length([arcintv{:}]));
timepointscharacterization=repmat(NaN,1,length([arcintv{:}]));
connectingtimepointsindex=zeros(1,numberofstages-1);
arc2part=repmat(NaN,1,length([arcintv{:}]));
partstructure=zeros(1,numberofstages);
totalarcarg=repmat(NaN,1,length([arcintv{:}]));
countertimepoints=0;
arccounter=0;
for ii=1:numberofstages
    countertimepoints0=countertimepoints+1;
    if ii==1
        countertimepoints=countertimepoints+length(arcintv{ii});
        timepoints(countertimepoints0:countertimepoints)=arcintv{ii};
        timepointscharacterization(countertimepoints0:countertimepoints)=1;
        timepointscharacterization(1)=0;
        if OCMATFTE.optimalswitchingtime(ii)
            timepointscharacterization(countertimepoints)=-1;
        else
            timepointscharacterization(countertimepoints)=0;
        end
    else
        countertimepoints=countertimepoints+length(arcintv{ii})-1;
        timepoints(countertimepoints0:countertimepoints)=arcintv{ii}(2:end);
        timepointscharacterization(countertimepoints0:countertimepoints)=1;
        if OCMATFTE.optimalswitchingtime(ii)
            timepointscharacterization(countertimepoints)=-1;
        else
            timepointscharacterization(countertimepoints)=0;
        end
    end
    if ii<numberofstages
        connectingtimepointsindex(ii)=countertimepoints;
    end
    arcarg=arcargument(ocTrjMP(ii));
    arccounter0=arccounter+1;
    for jj=1:length(arcarg)
        arccounter=arccounter+1;
        arc2part(arccounter)=ii;
    end
    partstructure(ii)=arccounter;
    totalarcarg(arccounter0:arccounter)=arcarg;
end
timepoints(isnan(timepoints))=[];
timepointscharacterization(isnan(timepointscharacterization))=[];
arc2part(isnan(arc2part))=[];
totalarcarg(isnan(totalarcarg))=[];

OCMATFTE.edge=[totalarcarg(1:end-1);totalarcarg(2:end)];
OCMATFTE.partstructure=partstructure;
OCMATFTE.initialtimepoints=timepoints;
OCMATFTE.timepointscharacterization=timepointscharacterization;
OCMATFTE.arc2part=arc2part;
OCMATFTE.partposition=[1,partstructure(1:end-1)+1;partstructure];
if ~isempty(connectiontimeindex)
    OCMATFTE.connectiontimeindex=connectiontimeindex;
else
    OCMATFTE.connectiontimeindex=connectingtimepointsindex;
end
switch continuationtype
    case 'initialstate'
        OCMATFTE.continuationvector=continuationtarget(:)-OCMATFTE.initialstate;
        continuationparameter=0;
    case 'endstate'
        OCMATFTE.continuationvector=continuationtarget(:)-OCMATFTE.endstate;
        continuationparameter=0;

    case {'switchingtime','time'}
        OCMATFTE.timepointscharacterization(continuationindex)=2;
        OCMATFTE.initialtime4continuation=timepoints(continuationindex);
        OCMATFTE.continuationindex=continuationindex;
        OCMATFTE.continuationindex2vector=1:length(continuationindex);
        OCMATFTE.continuationvector=continuationtarget(:).'-OCMATFTE.initialtime4continuation(:).';
        continuationparameter=0;
    case 'parameter'
        
        for ii=1:numberofstages
            OCMATFTE.initialparameter{ii}=OCMATFTE.parametervalue{ii}(continuationindex{ii});
            OCMATFTE.initialparameter{ii}=OCMATFTE.initialparameter{ii};
            if ~isempty(continuationtarget{ii})
                OCMATFTE.continuationvector{ii}=continuationtarget{ii}-OCMATFTE.initialparameter{ii};
            else
                OCMATFTE.continuationvector{ii}=[];
            end
            if isempty(continuationindex{ii})
                OCMATFTE.continuationvector{ii}=0;
            end
        end
        %         if iscell(continuationindex)
        %         elseif length(parindex{1})==1
        %             for ii=1:numberofstages
        %                 OCMATFTE.initialparameter(ii)=OCMATFTE.parametervalue{ii}(parindex{ii});
        %             end
        %             OCMATFTE.continuationvector=continuationtarget(:)-OCMATFTE.initialparameter(:);
        %             OCMATFTE.initialparameter=OCMATFTE.initialparameter(:);
        %         else
        %             OCMATFTE.initialparameter=OCMATFTE.parametervalue{1}(parindex{1});
        %             OCMATFTE.continuationvector=continuationtarget(:)-OCMATFTE.initialparameter(:);
        %             OCMATFTE.continuationvector=repmat(OCMATFTE.continuationvector(:).',OCMATFTE.numberofstages,1);
        %             OCMATFTE.initialparameter=repmat(OCMATFTE.initialparameter(:).',OCMATFTE.numberofstages,1);
        %         end
        continuationparameter=0;
end

% add free and optimized timepoints to free parameter values
OCMATFTE.timepointsfreeindex=find(OCMATFTE.timepointscharacterization==1);
OCMATFTE.timepointsfree4parcoordinate=length(sol.parameters)+(1:length(OCMATFTE.timepointsfreeindex));
sol.parameters=[sol.parameters timepoints(OCMATFTE.timepointscharacterization==1)];
OCMATFTE.timepointsoptimizedindex=find(OCMATFTE.timepointscharacterization==-1);
OCMATFTE.timepointsoptimize4parcoordinate=length(sol.parameters)+(1:length(OCMATFTE.timepointsoptimizedindex));
sol.parameters=[sol.parameters timepoints(OCMATFTE.timepointscharacterization==-1)];
OCMATFTE.continuationcoordinate=length(sol.parameters)+1;
sol.parameters=[sol.parameters continuationparameter];

OCMATFTE.continuationcoordinate=length(sol.parameters);
OCMATFTE.continuationtarget=continuationtarget;
OCMATFTE.continuationindex=continuationindex;

OCMATFTE.statecoordinate=scoord;
OCMATFTE.fixinitialstatecoord=fixinitialstate;
OCMATFTE.fixendstatecoord=fixendstate;
pathname=OCMATFTE.datapath{1}();
[resultfile,globalvarfile]=OCMATFTE.saveintermediatefiles{1}();
OCMATFTE.basicresultfilename=fullocmatfile(pathname,resultfile);
OCMATFTE.basicglobalvarfilename=fullocmatfile(pathname,globalvarfile);

OCMATCONT.HE.numinitialcondition=[];

OCMATCONT.codimension=1;

function sol=generatesolstruct(ocMultiPath)
global OCMATFTE


nummult=numberofparts(ocMultiPath);
sol.x=independentvar(ocMultiPath(1));
sol.y=dependentvar(ocMultiPath(1));
sol.arcarg=arcargument(ocMultiPath(1));
sol.arcinterval=arcinterval(ocMultiPath(1));
sol.parameters=parameters(ocMultiPath(1));
x0=initialtime(ocMultiPath(1));
for ii=2:nummult
    sol.x=[sol.x independentvar(ocMultiPath(ii))-x0+sol.x(end)];

    sol.y=[sol.y dependentvar(ocMultiPath(ii))];
    sol.arcarg=[sol.arcarg arcargument(ocMultiPath(ii))];
    actarcinterval=arcinterval(ocMultiPath(ii));
    sol.arcinterval=[sol.arcinterval actarcinterval(2:end)];
end
sol.parameters=[];
sol.x0=x0;
sol.solver='';
sol.parameters=[];
sol.solverinfo.tangent=[];
sol.solverinfo.coeff=[];
sol.solverinfo.tmesh=[];
