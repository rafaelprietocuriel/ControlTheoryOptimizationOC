function ocStruct=ppdemodelderivenecessaryconditions(ocStruct,symkernel,opt)
%
% derive the necessary optimality conditions from the initialization file
if nargin==2
    opt=defaultocoptions;
end
JacobianCalc=getocoptions(opt,'INIT','Jacobian');
HessianCalc=getocoptions(opt,'INIT','Hessian');
if ~isstruct(ocStruct)
    ocStruct=loadmodeldata(modelname(ocStruct));
end
% general properties
arcidentifier=retrieveppdemodelinformation(ocStruct,'arcidentifier');
inequalitycontrolconstraintnum=retrieveppdemodelinformation(ocStruct,'inequalitycontrolconstraintnum');
arcnum=retrieveppdemodelinformation(ocStruct,'arcnum');
statenum=retrieveppdemodelinformation(ocStruct,'statenum');

%
% The Pontryagin function and its derivatives
%
ocStruct.pontryaginfunction.identifier='PF';
P=retrieveppdemodelinformation(ocStruct,'pontryaginfunction');
Px=retrieveppdemodelinformation(ocStruct,'pontryaginfunctionDx',[],symkernel);
Pu=retrieveppdemodelinformation(ocStruct,'pontryaginfunctionDu',[],symkernel);
Px2=retrieveppdemodelinformation(ocStruct,'pontryaginfunctionDx2',[],symkernel);
Pu2=retrieveppdemodelinformation(ocStruct,'pontryaginfunctionDu2',[],symkernel);
Pux=retrieveppdemodelinformation(ocStruct,'pontryaginfunctionDuDx',[],symkernel);
ocStruct.pontryaginfunction.term=P.value;
ocStruct.pontryaginfunction.arcidentifier=arcidentifier.value;
ocStruct.pontryaginfunction.derivative.Dx.term=Px.value;
ocStruct.pontryaginfunction.derivative.Dx.identifier='DPFDx';
ocStruct.pontryaginfunction.derivative.Du.term=Pu.value;
ocStruct.pontryaginfunction.derivative.Du.identifier='DPFDu';
ocStruct.pontryaginfunction.derivative.Dx2.term=Px2.value;
ocStruct.pontryaginfunction.derivative.Dx2.identifier='DPFDx2';
ocStruct.pontryaginfunction.derivative.Du2.term=Pu2.value;
ocStruct.pontryaginfunction.derivative.Du2.identifier='DPFDu2';
ocStruct.pontryaginfunction.derivative.DuDx.term=Pux.value;
ocStruct.pontryaginfunction.derivative.DuDx.identifier='DPFDuDx';

%
% The adjoint system
%
% determine possible algebraic equations
dldt=retrieveppdemodelinformation(ocStruct,'adjointsystem',arcidentifier.value,symkernel);
statedynamicstype=retrieveppdemodelinformation(ocStruct,'statedynamicstype');
adjointsystemconvection=retrieveppdemodelinformation(ocStruct,'adjointsystemconvection');
adjointsystemdiffusion=retrieveppdemodelinformation(ocStruct,'adjointsystemdiffusion');
adjointsystemboundarycondition=retrieveppdemodelinformation(ocStruct,'adjointsystemboundarycondition',arcidentifier.value,symkernel);
for ii=1:statenum.value
    if ~isempty(adjointsystemconvection.value)
        ocStruct.foc.adjointsystem.dynamics.(arcidentifier2field(arcidentifier.value)).costate(ii).convection=adjointsystemconvection.value{ii};
    end
    if ~isempty(adjointsystemdiffusion.value)
        ocStruct.foc.adjointsystem.dynamics.(arcidentifier2field(arcidentifier.value)).costate(ii).diffusion=adjointsystemdiffusion.value{ii};
    end
    ocStruct.foc.adjointsystem.dynamics.(arcidentifier2field(arcidentifier.value)).costate(ii).term=dldt.value{ii};
    ocStruct.foc.adjointsystem.dynamics.(arcidentifier2field(arcidentifier.value)).costate(ii).type=statedynamicstype.value{ii};
    ocStruct.foc.adjointsystem.dynamics.(arcidentifier2field(arcidentifier.value)).costate(ii).boundarycondition=adjointsystemboundarycondition.value{ii};
end
%
%
% Hamiltonian maximizing condition
%
% solve the corresponding first order necessary conditions and let user
% decide which solutions should be used if multiple solutions occur
%
% test condition for the locally uniqueness of the solution for Hu=0.
% Therefore the models standard parameter values are used.
parametervalue=retrievemodelinformation(ocStruct,'parametervalue');
parametername=retrievemodelinformation(ocStruct,'parametername');
Hu2=[];
for ii=1:length(Pu2.value)
    Hu2=[Hu2;mystr2sym(['[' Pu2.value{ii} ']'])];
end
parval=num2cell(parametervalue.value);
Hu2=subs(Hu2,mystr2sym(parametername.value.'),parval);
eigval=eig(Hu2);
Hu2flag=1;
    if any(eigval==0)
    Hu2flag=0;
    fprintf('\n')
    flag=1;
    while flag
        answer=input('d^2 H/du^2 is singular for the standard parameter values. Interrupt y/n (y) : ','s');
        if isempty(answer)
            answer='y';
        end
        if strcmpi(answer,'y')
            ocmaterror('Control variables may not appear nonlinearly.\n')
        elseif strcmpi(answer,'n')
            flag=0;
        end
    end
end
arcwitheequalpontryaginfunction=retrieveppdemodelinformation(ocStruct,'arcwitheequalpontryaginfunction');
maximizingvariable=retrieveppdemodelinformation(ocStruct,'maximizingvariable');
solutionindex=retrieveppdemodelinformation(ocStruct,'solutionindex');
for ii=1:numel(arcwitheequalpontryaginfunction.value)
    arcid=arcwitheequalpontryaginfunction.value{ii}{1};
    if Hu2flag
        sol=retrieveppdemodelinformation(ocStruct,'maximizingregularnonlinearsolution',arcid,symkernel);
    end
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
    if strcmpi(JacobianCalc,'explicit')
        statename=retrieveppdemodelinformation(ocStruct,'statename',arcidentifier.value{ii});
        costatename=retrieveppdemodelinformation(ocStruct,'costatename',arcidentifier.value{ii});
        dxdt=retrieveppdemodelinformation(ocStruct,'specificstatedynamics',arcidentifier.value{ii},symkernel);
        dldt=retrieveppdemodelinformation(ocStruct,'specificadjointsystem',arcidentifier.value{ii},symkernel);
        dXdt=cell2vectorstring(term2optterm_ppdemodel(ocStruct,[dxdt.value dldt.value],arcidentifier.value{ii},opt));
        variable=cell2vectorstring([statename.value costatename.value]);
        ocStruct.foc.canonicalsystem.derivative.(arc).DX.term=string2cell(removematrixstring(ocmatjacobian(dXdt,variable,symkernel,opt)),'matrix');

        parametername=retrieveppdemodelinformation(ocStruct,'parametername',arcidentifier.value{ii});
        variable=cell2vectorstring(parametername.value);
        ocStruct.foc.canonicalsystem.derivative.(arc).DP.term=string2cell(removematrixstring(ocmatjacobian(dXdt,variable,symkernel,opt)),'matrix');
    elseif strcmpi(JacobianCalc,'numerical')
        ocStruct.foc.canonicalsystem.derivative.(arc).DX.term='';
        ocStruct.foc.canonicalsystem.derivative.(arc).DP.term='';
    else
        ocStruct.foc.canonicalsystem.derivative.(arc).DX.term='';
        ocStruct.foc.canonicalsystem.derivative.(arc).DP.term='';
    end
end

%
% Hessian of the canonical system
%
for ii=1:arcnum.value
    if strcmpi(HessianCalc,'explicit')
        arc=arcidentifier2field(arcidentifier.value{ii});
        H=retrieveppdemodelinformation(ocStruct,'canonicalsystemhessian',arcidentifier.value{ii},symkernel);
        ocStruct.foc.canonicalsystem.derivative.(arc).DX2.term=H.value;
        Jxp=retrieveppdemodelinformation(ocStruct,'canonicalsystemtotalhessian',arcidentifier.value{ii},symkernel);
        ocStruct.foc.canonicalsystem.derivative.(arc).DXP.term=Jxp.value;
    elseif strcmpi(HessianCalc,'numerical')
        ocStruct.foc.canonicalsystem.derivative.(arc).DX2.term='';
        ocStruct.foc.canonicalsystem.derivative.(arc).DXP.term='';
    else
        ocStruct.foc.canonicalsystem.derivative.(arc).DX2.term='';
        ocStruct.foc.canonicalsystem.derivative.(arc).DXP.term='';
    end
end
