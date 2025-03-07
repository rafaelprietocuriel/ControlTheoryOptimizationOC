function ocStruct=standarddiffmodelderivenecessaryconditions(ocStruct,symkernel,opt)
%
% derive the necessary optimality conditions from the initialization file

% general properties
arcidentifier=retrievediffmodelinformation(ocStruct,'arcidentifier');
inequalitycontrolconstraintnum=retrievediffmodelinformation(ocStruct,'inequalitycontrolconstraintnum');
inequalitystateconstraintnum=retrievediffmodelinformation(ocStruct,'inequalitystateconstraintnum');
arcnum=retrievediffmodelinformation(ocStruct,'arcnum');

%
% The Pontryagin function and its derivatives
%
ocStruct.pontryaginfunction.identifier='PF';
P=retrievediffmodelinformation(ocStruct,'pontryaginfunction');
Pxt=retrievediffmodelinformation(ocStruct,'pontryaginfunctionDxt',[],symkernel);
Putp1=retrievediffmodelinformation(ocStruct,'pontryaginfunctionDutp1',[],symkernel);
Pxt2=retrievediffmodelinformation(ocStruct,'pontryaginfunctionDxt2',[],symkernel);
Putp12=retrievediffmodelinformation(ocStruct,'pontryaginfunctionDutp12',[],symkernel);
Putp1xt=retrievediffmodelinformation(ocStruct,'pontryaginfunctionDutp1Dxt',[],symkernel);
ocStruct.pontryaginfunction.term=P.value;
ocStruct.pontryaginfunction.arcidentifier=arcidentifier.value;
ocStruct.pontryaginfunction.derivative.Dxt.term=Pxt.value;
ocStruct.pontryaginfunction.derivative.Dxt.identifier='DPFDxt';
ocStruct.pontryaginfunction.derivative.Dutp1.term=Putp1.value;
ocStruct.pontryaginfunction.derivative.Dutp1.identifier='DPFDutp1';
ocStruct.pontryaginfunction.derivative.Dxt2.term=Pxt2.value;
ocStruct.pontryaginfunction.derivative.Dxt2.identifier='DPFDxt2';
ocStruct.pontryaginfunction.derivative.Dutp12.term=Putp12.value;
ocStruct.pontryaginfunction.derivative.Dutp12.identifier='DPFDutp12';
ocStruct.pontryaginfunction.derivative.Dutp1Dxt.term=Putp1xt.value;
ocStruct.pontryaginfunction.derivative.Dutp1Dxt.identifier='DPFDutp1Dxt';
if implicitcontrols(ocStruct)
    Pulmmc=retrievediffmodelinformation(ocStruct,'pontryaginfunctionDuDlmmc',[],symkernel);
    ocStruct.pontryaginfunction.derivative.DuDlmmc.term=Pulmmc.value;
    ocStruct.pontryaginfunction.derivative.DuDlmmc.identifier='DPFDuDlmmc';
    PuX=retrievediffmodelinformation(ocStruct,'pontryaginfunctionDuDX',[],symkernel);
    ocStruct.pontryaginfunction.derivative.DuDX.term=PuX.value;
    ocStruct.pontryaginfunction.derivative.DuDX.identifier='DPFDuDX';
end

%
% The adjoint system
%
% determine possible algebraic equations
dldt=retrievediffmodelinformation(ocStruct,'adjointsystem',arcidentifier.value,symkernel);
arcid=arcidentifier.value;
for ii=arcnum.value:-1:1
    ae=retrievediffmodelinformation(ocStruct,'algebraicequation',arcid{ii},symkernel);
    ocStruct.foc.adjointsystem.algebraicequation.(arcidentifier2field(arcid{ii})).term=ae.value;
end
ocStruct.foc.adjointsystem.dynamics.(arcidentifier2field(arcidentifier.value)).diff.term=dldt.value;
ocStruct.foc.adjointsystem.dynamics.(arcidentifier2field(arcidentifier.value)).diff.maptype='implicit';
ocStruct.foc.adjointsystem.algebraicequation=orderfields(ocStruct.foc.adjointsystem.algebraicequation);

%
%
% Hamiltonian maximizing condition
%
% solve the corresponding first order necessary conditions and let user
% decide which solutions should be used if multiple solutions occur
%
arcwitheequalpontryaginfunction=retrievediffmodelinformation(ocStruct,'arcwitheequalpontryaginfunction');
maximizingvariable=retrievediffmodelinformation(ocStruct,'maximizingvariable');
solutionindex=retrievediffmodelinformation(ocStruct,'solutionindex');
for ii=1:numel(arcwitheequalpontryaginfunction.value)
    arcid=arcwitheequalpontryaginfunction.value{ii}{1};
    sol=retrievediffmodelinformation(ocStruct,'maximizingsolution',arcid,symkernel);
    if numel(sol.value)>1
        solutionindices=cell2vectorstring(solutionindex.value{ii});
        % multiple solutions exist
        fprintf('Arc identifiers : %s\n\n',arcidentifier2field(arcwitheequalpontryaginfunction.value{ii}))
        for jj=1:numel(sol.value)
            fprintf('Solution Number %d\n',jj)
            for kk=1:numel(maximizingvariable.value)
                fprintf('\t%s = %s\n',maximizingvariable.value{kk},sol.value(jj).(maximizingvariable.value{kk}))
            end
            fprintf('\n\n')
        end
        fprintf('%d solutions can be kept?\nSolutions chosen by default are %s\n ',numel(solutionindex.value{ii}),solutionindices)
        while 1
            accept=lower(input('Accept: (y)/n : ','s'));
            if isempty(accept)
                accept='y';
            end
            if strcmp(accept,'y')
                sol.value=eval(['sol.value(' solutionindices ')']);
                break
            elseif strcmp(accept,'n')
                break
            elseif regexp(accept,'q')
                break
            end
        end
        counter=0;
        while 1 && ~strcmp(accept,'y')
            counter=counter+1;
            fprintf('\n\n')
            answer=input('Input solution indices in Matlab notation : ','s');
            solutionindices=eval(['[' answer ']']);
            try
                if numel(solutionindices)~=numel(solutionindex.value{ii})
                    fprintf('Wrong number of solutions. Try again.\n')
                elseif ~all(ismember(solutionindices,1:numel(sol.value)))
                    fprintf('Solution indices are out of solution range.\n')
                else
                    sol.value=sol.value(solutionindices);
                    break
                end
            catch
                fprintf('\nWrong format %s. Try again.\n',answer)
            end
            if counter>2
                fprintf('\nGive up after %d try.\n',counter)
            end
        end
    end
    for jj=1:numel(arcwitheequalpontryaginfunction.value{ii})
        ocStruct.foc.value.control.(arcidentifier2field(arcwitheequalpontryaginfunction.value{ii}{jj})).term=sol.value(jj).control;
        if inequalitycontrolconstraintnum.value
            ocStruct.foc.value.lagrangemultcc.(arcidentifier2field(arcwitheequalpontryaginfunction.value{ii}{jj})).term=sol.value(jj).lagrangemultcc;
        end
    end
end
% sort arc identifiers
arcargstring=char(fieldnames(ocStruct.foc.value.control));
[dum idx]=sort(str2num(arcargstring(:,4:end)));
ocStruct.foc.value.control=orderfields(ocStruct.foc.value.control,idx);
if inequalitycontrolconstraintnum.value
    ocStruct.foc.value.lagrangemultcc=orderfields(ocStruct.foc.value.lagrangemultcc,idx);
end

%
% Jacobian of the canonical system
%
for ii=1:arcnum.value
    arc=arcidentifier2field(arcidentifier.value{ii});
    J=retrievediffmodelinformation(ocStruct,'canonicalsystemjacobian',arcidentifier.value{ii},symkernel);
    ocStruct.foc.canonicalsystem.derivative.(arc).DX.term=J.value;
    Jp=retrievediffmodelinformation(ocStruct,'canonicalsystemparameterjacobian',arcidentifier.value{ii},symkernel);
    ocStruct.foc.canonicalsystem.derivative.(arc).DP.term=Jp.value;
end

%
% Hessian of the canonical system
%
% for ii=1:arcnum.value
%     arc=arcidentifier2field(arcidentifier.value{ii});
%     H=retrievediffmodelinformation(ocStruct,'canonicalsystemhessian',arcidentifier.value{ii},symkernel);
%     ocStruct.foc.canonicalsystem.derivative.(arc).DX2.term=H.value;
%     Jxp=retrievediffmodelinformation(ocStruct,'canonicalsystemtotalhessian',arcidentifier.value{ii},symkernel);
%     ocStruct.foc.canonicalsystem.derivative.(arc).DXP.term=Jxp.value;
% end