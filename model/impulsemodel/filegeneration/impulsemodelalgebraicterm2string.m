function algebraic=impulsemodelalgebraicterm2string(ocStruct,class,vectflag,varargin)
%
% ALGEBRAICTERM2STRING process the algebraic terms of the structure
% ocStruct describing the model into strings which are printed into the
% model files


algebraic.term='';
algebraic.type='';
algebraic.info=[];
arcfield='';
classtype='';

if isempty(ocStruct)
    return
end

if nargin==2
    vectflag=0;
end

if nargin>=4
    arcfield=varargin{1};
end

if nargin>=5
    classtype=varargin{2};
end

switch class
    case 'controlvalue'
        optimalvalue=retrieveimpulsemodelinformation(ocStruct,'optimalvalue',field2arcidentifier(arcfield));
        controlnum=retrieveimpulsemodelinformation(ocStruct,'controlnum');
        controlname=retrieveimpulsemodelinformation(ocStruct,'controlname');
        if controlnum.value
            for ii=1:controlnum.value
                if ii==1
                    if controlnum.value>1
                        algebraic.term{ii}=['[' optimalvalue.value.(controlname.value{ii}) '; ...'];
                    else
                        algebraic.term{ii}=['[' optimalvalue.value.(controlname.value{ii})];
                    end
                elseif ii<controlnum.value
                    algebraic.term{ii}=[optimalvalue.value.(controlname.value{ii}) '; ...'];
                else
                    algebraic.term{ii}=[optimalvalue.value.(controlname.value{ii})];
                end
            end
            algebraic.term{ii}=[algebraic.term{ii} ']'];
        else
            algebraic.term='';
        end

    case 'impulsecontrolvalue'
        ioptimalvalue=retrieveimpulsemodelinformation(ocStruct,'impulseoptimalvalue',field2arcidentifier(arcfield));
        icontrolnum=retrieveimpulsemodelinformation(ocStruct,'icontrolnum');
        icontrolname=retrieveimpulsemodelinformation(ocStruct,'icontrolname');
        for ii=1:icontrolnum.value
            if ii==1
                if icontrolnum.value>1
                    algebraic.term{ii}=['[' ioptimalvalue.value.(icontrolname.value{ii}) '; ...'];
                else
                    algebraic.term{ii}=[ioptimalvalue.value.(icontrolname.value{ii})];
                end
            elseif ii<icontrolnum.value
                algebraic.term{ii}=[ioptimalvalue.value.(icontrolname.value{ii}) '; ...'];
            else
                algebraic.term{ii}=[ioptimalvalue.value.(icontrolname.value{ii})];
            end
        end
        if icontrolnum.value>1
            algebraic.term{ii}=[algebraic.term{ii} ']'];
        end
    case 'symbolicimpulsecontrolvalue'
        ioptimalvalue=retrieveimpulsemodelinformation(ocStruct,'impulseoptimalvalue',field2arcidentifier(arcfield));
        icontrolnum=retrieveimpulsemodelinformation(ocStruct,'icontrolnum');
        icontrolname=retrieveimpulsemodelinformation(ocStruct,'icontrolname');
        for ii=1:icontrolnum.value
            if ii==1
                if icontrolnum.value>1
                    algebraic.term{ii}=['sym([''[' ioptimalvalue.value.(icontrolname.value{ii}) ','' ...'];
                else
                    algebraic.term{ii}=['sym([''[' ioptimalvalue.value.(icontrolname.value{ii}) ']''])'];
                end
            elseif ii<icontrolnum.value
                algebraic.term{ii}=['''' ioptimalvalue.value.(icontrolname.value{ii}) ','' ...'];
            else
                algebraic.term{ii}=['''' ioptimalvalue.value.(icontrolname.value{ii}) ']''])'];
            end
        end

    case 'optimalvalue4statecontrol'
        optimalvalue=retrieveimpulsemodelinformation(ocStruct,'optimalvalue4statecontrol',field2arcidentifier(arcfield));
        controlnum=retrieveimpulsemodelinformation(ocStruct,'controlnum');
        controlname=retrieveimpulsemodelinformation(ocStruct,'controlname');
        for ii=1:controlnum.value
            if ii==1
                if controlnum.value>1
                    algebraic.term{ii}=['[' optimalvalue.value.(controlname.value{ii}) '; ...'];
                else
                    algebraic.term{ii}=['[' optimalvalue.value.(controlname.value{ii})];
                end
            elseif ii<controlnum.value
                algebraic.term{ii}=[optimalvalue.value.(controlname.value{ii}) '; ...'];
            else
                algebraic.term{ii}=[optimalvalue.value.(controlname.value{ii})];
            end
        end
        algebraic.term{ii}=[algebraic.term{ii} ']'];

    case 'optimalcostatevalue'
        optimalcostatevalue=retrieveimpulsemodelinformation(ocStruct,'optimalcostatevalue',field2arcidentifier(arcfield),getsymkernel);
        costatenum=retrieveimpulsemodelinformation(ocStruct,'costatenum');
        costatename=retrieveimpulsemodelinformation(ocStruct,'costatename');
        for ii=1:costatenum.value
            if ii==1
                if costatenum.value>1
                    algebraic.term{ii}=['[' optimalcostatevalue.value.(costatename.value{ii}) '; ...'];
                else
                    algebraic.term{ii}=optimalcostatevalue.value.(costatename.value{ii});
                end
            elseif ii<costatenum.value
                algebraic.term{ii}=[optimalcostatevalue.value.(costatename.value{ii}) '; ...'];
            else
                algebraic.term{ii}=[optimalcostatevalue.value.(costatename.value{ii})];
            end
        end
        if costatenum.value>1
            algebraic.term{ii}=[algebraic.term{ii} ']'];
        end

    case 'explicitstatevalue'
        explicitstatevalue=retrieveimpulsemodelinformation(ocStruct,'explicitstatevalue',field2arcidentifier(arcfield),getsymkernel);
        statenum=retrieveimpulsemodelinformation(ocStruct,'statenum');
        for ii=1:statenum.value
            if ii==1
                if statenum.value>1
                    algebraic.term{ii}=['[' explicitstatevalue.value{ii} '; ...'];
                else
                    algebraic.term{ii}=explicitstatevalue.value{ii};
                end
            elseif ii<statenum.value
                algebraic.term{ii}=[explicitstatevalue.value{ii} '; ...'];
            else
                algebraic.term{ii}=[explicitstatevalue.value{ii}];
            end
        end
        if statenum.value>1
            algebraic.term{ii}=[algebraic.term{ii} ']'];
        end

    case 'symbolicstatevalue'
        explicitstatevalue=retrieveimpulsemodelinformation(ocStruct,'explicitstatevalue',field2arcidentifier(arcfield),getsymkernel);
        statenum=retrieveimpulsemodelinformation(ocStruct,'statenum');
        for ii=1:statenum.value
            if ii==1
                if statenum.value>1
                    algebraic.term{ii}=['sym([''[' explicitstatevalue.value{ii} ','' ...'];
                else
                    algebraic.term{ii}=['sym([''[' explicitstatevalue.value{ii} ']''])'];
                end
            elseif ii<statenum.value
                algebraic.term{ii}=['''' explicitstatevalue.value{ii} ','' ...'];
            else
                algebraic.term{ii}=['''' explicitstatevalue.value{ii} ']''])'];
            end
        end

    case 'equilibriumequation'
        equilibriumequation=retrieveimpulsemodelinformation(ocStruct,'equilibriumequation',field2arcidentifier(arcfield),getsymkernel);
        equationnum=length(equilibriumequation.value);
        for ii=1:equationnum
            if ii==1
                if equationnum>1
                    algebraic.term{ii}=['sym([''[' equilibriumequation.value{ii} ','' ...'];
                else
                    algebraic.term{ii}=['sym([''[' equilibriumequation.value{ii} ']''])'];
                end
            elseif ii<equationnum
                algebraic.term{ii}=['''' equilibriumequation.value{ii} ','' ...'];
            else
                algebraic.term{ii}=['''' equilibriumequation.value{ii} ']''])'];
            end
        end

    case 'symboliccontrolvalue'
        optimalvalue=retrieveimpulsemodelinformation(ocStruct,'optimalvalue',field2arcidentifier(arcfield));
        controlnum=retrieveimpulsemodelinformation(ocStruct,'controlnum');
        controlname=retrieveimpulsemodelinformation(ocStruct,'controlname');
        for ii=1:controlnum.value
            if ii==1
                if controlnum.value>1
                    algebraic.term{ii}=['sym([''[' optimalvalue.value.(controlname.value{ii}) ','' ...'];
                else
                    algebraic.term{ii}=['sym([''[' optimalvalue.value.(controlname.value{ii}) ']''])'];
                end
            elseif ii<controlnum.value
                algebraic.term{ii}=['''' optimalvalue.value.(controlname.value{ii}) ','' ...'];
            else
                algebraic.term{ii}=['''' optimalvalue.value.(controlname.value{ii}) ']''])'];
            end
        end

    case 'symboliccontrolvalue4statecontrol'
        optimalvalue=retrieveimpulsemodelinformation(ocStruct,'optimalvalue4statecontrol',field2arcidentifier(arcfield));
        controlnum=retrieveimpulsemodelinformation(ocStruct,'controlnum');
        controlname=retrieveimpulsemodelinformation(ocStruct,'controlname');
        for ii=1:controlnum.value
            if ii==1
                if controlnum.value>1
                    algebraic.term{ii}=['sym([''[' optimalvalue.value.(controlname.value{ii}) ','' ...'];
                else
                    algebraic.term{ii}=['sym([''[' optimalvalue.value.(controlname.value{ii}) ']''])'];
                end
            elseif ii<controlnum.value
                algebraic.term{ii}=['''' optimalvalue.value.(controlname.value{ii}) ','' ...'];
            else
                algebraic.term{ii}=['''' optimalvalue.value.(controlname.value{ii}) ']''])'];
            end
        end
    case 'symboliccostatevalue'
        optimalvalue=retrieveimpulsemodelinformation(ocStruct,'optimalcostatevalue',field2arcidentifier(arcfield),getsymkernel);
        costatenum=retrieveimpulsemodelinformation(ocStruct,'costatenum');
        costatename=retrieveimpulsemodelinformation(ocStruct,'costatename');
        for ii=1:costatenum.value
            if ii==1
                if costatenum.value>1
                    algebraic.term{ii}=['sym([''[' optimalvalue.value.(costatename.value{ii}) ','' ...'];
                else
                    algebraic.term{ii}=['sym([''[' optimalvalue.value.(costatename.value{ii}) ']''])'];
                end
            elseif ii<costatenum.value
                algebraic.term{ii}=['''' optimalvalue.value.(costatename.value{ii}) ','' ...'];
            else
                algebraic.term{ii}=['''' optimalvalue.value.(costatename.value{ii}) ']''])'];
            end
        end

    case 'lagrangemultcc'
        optimalvalue=retrieveimpulsemodelinformation(ocStruct,'optimalvalue',field2arcidentifier(arcfield));
        inequalitycontrolconstraintnum=retrieveimpulsemodelinformation(ocStruct,'inequalitycontrolconstraintnum');
        lagrangemultipliercontrolname=retrieveimpulsemodelinformation(ocStruct,'lagrangemultipliercontrolname');
        for ii=1:inequalitycontrolconstraintnum.value
            if ii==1
                if inequalitycontrolconstraintnum.value>1
                    algebraic.term{ii}=['[' optimalvalue.value.(lagrangemultipliercontrolname.value{ii}) '; ...'];
                else
                    algebraic.term{ii}=['[' optimalvalue.value.(lagrangemultipliercontrolname.value{ii})];
                end
            elseif ii<inequalitycontrolconstraintnum.value
                algebraic.term{ii}=[optimalvalue.value.(lagrangemultipliercontrolname.value{ii}) '; ...'];
            else
                algebraic.term{ii}=optimalvalue.value.(lagrangemultipliercontrolname.value{ii});
            end
        end
        algebraic.term{ii}=[algebraic.term{ii} ']'];

    case 'symboliclagrangemultcc'
        optimalvalue=retrieveimpulsemodelinformation(ocStruct,'optimalvalue',field2arcidentifier(arcfield));
        inequalitycontrolconstraintnum=retrieveimpulsemodelinformation(ocStruct,'inequalitycontrolconstraintnum');
        lagrangemultipliercontrolname=retrieveimpulsemodelinformation(ocStruct,'lagrangemultipliercontrolname');
        for ii=1:inequalitycontrolconstraintnum.value
            if ii==1
                if inequalitycontrolconstraintnum.value>1
                    algebraic.term{ii}=['sym([''[' optimalvalue.value.(lagrangemultipliercontrolname.value{ii}) ','' ...'];
                else
                    algebraic.term{ii}=['sym([''[' optimalvalue.value.(lagrangemultipliercontrolname.value{ii}) ']''])'];
                end

            elseif ii<inequalitycontrolconstraintnum.value
                algebraic.term{ii}=['''' optimalvalue.value.(lagrangemultipliercontrolname.value{ii}) ','' ...'];
            else
                algebraic.term{ii}=['''' optimalvalue.value.(lagrangemultipliercontrolname.value{ii}) ']''])'];
            end
        end

    case 'lagrangemultsc'
        optimalvalue=retrieveimpulsemodelinformation(ocStruct,'optimalvalue',field2arcidentifier(arcfield));
        inequalitystateconstraintnum=retrieveimpulsemodelinformation(ocStruct,'inequalitystateconstraintnum');
        lagrangemultiplierstatename=retrieveimpulsemodelinformation(ocStruct,'lagrangemultiplierstatename');
        for ii=1:inequalitystateconstraintnum.value
            if ii==1
                if inequalitystateconstraintnum.value>1
                    algebraic.term{ii}=['[' optimalvalue.value.(lagrangemultiplierstatename.value{ii}) '; ...'];
                else
                    algebraic.term{ii}=['[' optimalvalue.value.(lagrangemultiplierstatename.value{ii})];
                end
            elseif ii<inequalitystateconstraintnum.value
                algebraic.term{ii}=[optimalvalue.value.(lagrangemultiplierstatename.value{ii}) '; ...'];
            else
                algebraic.term{ii}=optimalvalue.value.(lagrangemultiplierstatename.value{ii});
            end
        end
        algebraic.term{ii}=[algebraic.term{ii} ']'];

    case 'symboliclagrangemultsc'
        optimalvalue=retrieveimpulsemodelinformation(ocStruct,'optimalvalue',field2arcidentifier(arcfield));
        inequalitystateconstraintnum=retrieveimpulsemodelinformation(ocStruct,'inequalitystateconstraintnum');
        lagrangemultiplierstatename=retrieveimpulsemodelinformation(ocStruct,'lagrangemultiplierstatename');
        for ii=1:inequalitystateconstraintnum.value
            if ii==1
                if inequalitystateconstraintnum.value>1
                    algebraic.term{ii}=['sym([''[' optimalvalue.value.(lagrangemultiplierstatename.value{ii}) ','' ...'];
                else
                    algebraic.term{ii}=['sym([''[' optimalvalue.value.(lagrangemultiplierstatename.value{ii}) ']''])'];
                end

            elseif ii<inequalitystateconstraintnum.value
                algebraic.term{ii}=['''' optimalvalue.value.(lagrangemultiplierstatename.value{ii}) ','' ...'];
            else
                algebraic.term{ii}=['''' optimalvalue.value.(lagrangemultiplierstatename.value{ii}) ']''])'];
            end
        end

    case 'objective'
        objectiveintegrand=retrieveimpulsemodelinformation(ocStruct,'objectiveintegrand');
        switch classtype
            case 'current'
                discountfactor=retrieveimpulsemodelinformation(ocStruct,'discountfactor');
                algebraic.term=[discountfactor.value '*(' objectiveintegrand.value ')'];
            case 'present'
                algebraic.term=objectiveintegrand.value;
        end

    case 'salvagevalue'
        salvagevalue=retrieveimpulsemodelinformation(ocStruct,'salvagevalue');
        %endtimediscountfactor=retrieveimpulsemodelinformation(ocStruct,'endtimediscountfactor');
        algebraic.term=salvagevalue.value;

    case 'discountedsalvagevalue'
        salvagevalue=retrieveimpulsemodelinformation(ocStruct,'salvagevalue');
        endtimediscountfactor=retrieveimpulsemodelinformation(ocStruct,'endtimediscountfactor');
        algebraic.term=[endtimediscountfactor.value '*(' salvagevalue.value ')'];

    case 'objectivesummand'
        objectivesummand=retrieveimpulsemodelinformation(ocStruct,'objectivesummand');
        algebraic.term=objectivesummand.value;

    case 'discountedobjectivesummand'
        objectivesummand=retrieveimpulsemodelinformation(ocStruct,'discountedobjectivesummand');
        algebraic.term=objectivesummand.value;

    case 'controlconstraint'
        switch classtype
            case 'inequality'
                inequalitycontrolconstraint=retrieveimpulsemodelinformation(ocStruct,'inequalitycontrolconstraint');
                algebraic.term=inequalitycontrolconstraint.value;
            otherwise
        end

    case 'constraint'
        inequalitycontrolconstraint=retrieveimpulsemodelinformation(ocStruct,'inequalitycontrolconstraint');
        inequalitystateconstraint=retrieveimpulsemodelinformation(ocStruct,'inequalitystateconstraint');
        algebraic.term=[inequalitycontrolconstraint.value inequalitystateconstraint.value];

    case 'canonicalsystem'
        canonicalsystem=retrieveimpulsemodelinformation(ocStruct,'canonicalsystem','','');
        for ii=1:numel(canonicalsystem.value)
            if ii==1
                algebraic.term{ii}=['[' canonicalsystem.value{ii} '; ...'];
            elseif ii<numel(canonicalsystem.value)
                algebraic.term{ii}=[canonicalsystem.value{ii} '; ...'];
            else
                algebraic.term{ii}=[canonicalsystem.value{ii} ']'];
            end
        end

    case 'exogenousfunctionterm'
        exogenousfunctionterm=retrieveimpulsemodelinformation(ocStruct,'exogenousfunctionterm');
        exogenousfunctionnum=retrieveimpulsemodelinformation(ocStruct,'exogenousfunctionnum');
        for ii=1:exogenousfunctionnum.value
            if ii==1
                if exogenousfunctionnum.value>1
                    algebraic.term{ii}=['out=[' exogenousfunctionterm.value{ii} '; ...'];
                else
                    algebraic.term{ii}=['out=' exogenousfunctionterm.value{ii} ';'];
                end
            elseif ii<exogenousfunctionnum.value
                algebraic.term{ii}=[exogenousfunctionterm.value{ii} '; ...'];
            else
                algebraic.term{ii}=[exogenousfunctionterm.value{ii} ']'];
            end
        end

    case 'algebraicequation'
        arcidentifier=field2arcidentifier(arcfield);
        algebraicequation=retrieveimpulsemodelinformation(ocStruct,'algebraicequation',arcidentifier,getsymkernel);
        algebraicequationnum=retrieveimpulsemodelinformation(ocStruct,'algebraicequationnum',arcidentifier,'');
        for ii=1:numel(algebraicequation.value)
            if ii==1
                if algebraicequationnum.value>1
                    algebraic.term{ii}=['[' algebraicequation.value{ii} '; ...'];
                else
                    algebraic.term{ii}=algebraicequation.value{ii};
                end
            elseif ii<numel(algebraicequation.value)
                algebraic.term{ii}=['' algebraicequation.value{ii} '; ...'];
            else
                algebraic.term{ii}=['' algebraicequation.value{ii} ']'];
            end
        end

    case 'algebraicequationimplicit'
        arcidentifier=field2arcidentifier(arcfield);
        algebraicequationimplicit=retrieveimpulsemodelinformation(ocStruct,'algebraicequationimplicit',arcidentifier,getsymkernel);
        algebraicequationimplicitnum=numel(algebraicequationimplicit.value);
        for ii=1:algebraicequationimplicitnum
            if ii==1
                if algebraicequationimplicitnum>1
                    algebraic.term{ii}=['[' algebraicequationimplicit.value{ii} '; ...'];
                else
                    algebraic.term{ii}=algebraicequationimplicit.value{ii};
                end
            elseif ii<numel(algebraicequationimplicit.value)
                algebraic.term{ii}=['' algebraicequationimplicit.value{ii} '; ...'];
            else
                algebraic.term{ii}=['' algebraicequationimplicit.value{ii} ']'];
            end
        end

    case 'symboliccanonicalsystem'
        canonicalsystem=retrieveimpulsemodelinformation(ocStruct,'canonicalsystem','','');
        for ii=1:numel(canonicalsystem.value)
            if ii==1
                algebraic.term{ii}=['sym([''[' canonicalsystem.value{ii} ';'' ...'];
            elseif ii<numel(canonicalsystem.value)
                algebraic.term{ii}=['''' canonicalsystem.value{ii} ';'' ...'];
            else
                algebraic.term{ii}=['''' canonicalsystem.value{ii} ']''])'];
            end
        end

    case 'symbolicalgebraicequation'
        arcidentifier=field2arcidentifier(arcfield);
        algebraicequation=retrieveimpulsemodelinformation(ocStruct,'algebraicequation',arcidentifier,getsymkernel);
        algebraicequationnum=retrieveimpulsemodelinformation(ocStruct,'algebraicequationnum',arcidentifier,'');
        for ii=1:numel(algebraicequation.value)
            if ii==1
                if algebraicequationnum.value>1
                    algebraic.term{ii}=['sym([''[' algebraicequation.value{ii} ';'' ...'];
                else
                    algebraic.term{ii}=['sym(''' algebraicequation.value{ii} ''')'];
                end
            elseif ii<numel(algebraicequation.value)
                algebraic.term{ii}=['''' algebraicequation.value{ii} ';'' ...'];
            else
                algebraic.term{ii}=['''' algebraicequation.value{ii} ']''])'];
            end
        end

    case 'symboliccanonicalsystemjacobian'
        arcidentifier=field2arcidentifier(arcfield);
        canonicalsystemjacobian=retrieveimpulsemodelinformation(ocStruct,'canonicalsystemjacobian',arcidentifier,getsymkernel);
        numterm=numel(canonicalsystemjacobian.value);
        for ii=1:numterm
            if ii==1
                algebraic.term{ii}=['sym([''' canonicalsystemjacobian.value{ii} ';'' ...'];
            elseif ii<numterm
                algebraic.term{ii}=['''' canonicalsystemjacobian.value{ii} ';'' ...'];
            else
                algebraic.term{ii}=['''' canonicalsystemjacobian.value{ii} ''''];
            end
        end
        algebraic.term{ii}=[algebraic.term{ii} '])'];

    case 'symbolicstatecontrolcanonicalsystemjacobian'
        arcidentifier=field2arcidentifier(arcfield);
        canonicalsystemjacobian=retrieveimpulsemodelinformation(ocStruct,'canonicalsystemjacobianstatecontrol',arcidentifier,getsymkernel);
        numterm=numel(canonicalsystemjacobian.value);
        for ii=1:numterm
            if ii==1
                algebraic.term{ii}=['sym([''' canonicalsystemjacobian.value{ii} ';'' ...'];
            elseif ii<numterm
                algebraic.term{ii}=['''' canonicalsystemjacobian.value{ii} ';'' ...'];
            else
                algebraic.term{ii}=['''' canonicalsystemjacobian.value{ii} ''''];
            end
        end
        algebraic.term{ii}=[algebraic.term{ii} '])'];

    case 'symbolicstatecontrolcanonicalsystem'
        arcidentifier=field2arcidentifier(arcfield);
        specificcanonicalsysteminstatecontrol=retrieveimpulsemodelinformation(ocStruct,'specificcanonicalsysteminstatecontrol',arcidentifier,getsymkernel);
        numterm=numel(specificcanonicalsysteminstatecontrol.value);
        for ii=1:numterm
            if ii==1
                algebraic.term{ii}=['sym([''[' specificcanonicalsysteminstatecontrol.value{ii} ';'' ...'];
            elseif ii<numterm
                algebraic.term{ii}=['''' specificcanonicalsysteminstatecontrol.value{ii} ';'' ...'];
            else
                algebraic.term{ii}=['''' specificcanonicalsysteminstatecontrol.value{ii} ']'''];
            end
        end
        algebraic.term{ii}=[algebraic.term{ii} '])'];

    case 'specificcanonicalsystem'
        arcidentifier=field2arcidentifier(arcfield);
        specificcanonicalsystem=retrieveimpulsemodelinformation(ocStruct,'specificcanonicalsystem',arcidentifier,getsymkernel);
        numterm=numel(specificcanonicalsystem.value);
        for ii=1:numterm
            if ii==1
                algebraic.term{ii}=['[' specificcanonicalsystem.value{ii} '; ...'];
            elseif ii<numterm
                algebraic.term{ii}=[specificcanonicalsystem.value{ii} '; ...'];
            else
                algebraic.term{ii}=specificcanonicalsystem.value{ii};
            end
        end
        algebraic.term{ii}=[algebraic.term{ii} ']'];

    case 'specificpontryaginfunction'
        arcidentifier=field2arcidentifier(arcfield);
        specificpontryaginfunction=retrieveimpulsemodelinformation(ocStruct,'specificpontryaginfunction',arcidentifier,getsymkernel);

        algebraic.term=specificpontryaginfunction.value;

    case 'symboliccanonicalsystem'
        arcidentifier=field2arcidentifier(arcfield);
        specificcanonicalsystem=retrieveimpulsemodelinformation(ocStruct,'symboliccanonicalsystem',arcidentifier,getsymkernel);
        numterm=numel(specificcanonicalsystem.value);
        for ii=1:numterm
            if ii==1
                algebraic.term{ii}=['[' specificcanonicalsystem.value{ii} '; ...'];
            elseif ii<numterm
                algebraic.term{ii}=[specificcanonicalsystem.value{ii} '; ...'];
            else
                algebraic.term{ii}=specificcanonicalsystem.value{ii};
            end
        end
        algebraic.term{ii}=[algebraic.term{ii} ']'];

    case 'statecostatejacobian'
        arcidentifier=field2arcidentifier(arcfield);
        canonicalsystemjacobian=retrieveimpulsemodelinformation(ocStruct,'statecostatejacobian',arcidentifier,getsymkernel);
        numterm=numel(canonicalsystemjacobian.value);
        for ii=1:numterm
            if ii==1
                algebraic.term{ii}=['[' canonicalsystemjacobian.value{ii} '; ...'];
            elseif ii<numterm
                algebraic.term{ii}=[canonicalsystemjacobian.value{ii} '; ...'];
            else
                algebraic.term{ii}=canonicalsystemjacobian.value{ii};
            end
        end
        algebraic.term{ii}=[algebraic.term{ii} ']'];

    case 'canonicalsystemjacobian'
        arcidentifier=field2arcidentifier(arcfield);
        canonicalsystemjacobian=retrieveimpulsemodelinformation(ocStruct,'canonicalsystemjacobian',arcidentifier,getsymkernel);
        numterm=numel(canonicalsystemjacobian.value);
        for ii=1:numterm
            if ii==1
                algebraic.term{ii}=['[' canonicalsystemjacobian.value{ii} '; ...'];
            elseif ii<numterm
                algebraic.term{ii}=[canonicalsystemjacobian.value{ii} '; ...'];
            else
                algebraic.term{ii}=canonicalsystemjacobian.value{ii};
            end
        end
        algebraic.term{ii}=[algebraic.term{ii} ']'];

    case 'canonicalsystemjacobianstatecontrol'
        arcidentifier=field2arcidentifier(arcfield);
        canonicalsystemjacobianstatecontrol=retrieveimpulsemodelinformation(ocStruct,'canonicalsystemjacobianstatecontrol',arcidentifier,getsymkernel);
        numterm=numel(canonicalsystemjacobianstatecontrol.value);
        for ii=1:numterm
            if ii==1
                algebraic.term{ii}=['[' canonicalsystemjacobianstatecontrol.value{ii} '; ...'];
            elseif ii<numterm
                algebraic.term{ii}=[canonicalsystemjacobianstatecontrol.value{ii} '; ...'];
            else
                algebraic.term{ii}=canonicalsystemjacobianstatecontrol.value{ii};
            end
        end
        algebraic.term{ii}=[algebraic.term{ii} ']'];

    case 'optimalcontroldynamicsleftside'
        arcidentifier=field2arcidentifier(arcfield);
        optimalcontroldynamicsleftside=retrieveimpulsemodelinformation(ocStruct,'optimalcontroldynamicsleftside',arcidentifier,getsymkernel);
        implicitnonlinearcontrolindex=retrieveimpulsemodelinformation(ocStruct,'totalimplicitvariableindex');
        numterm=numel(optimalcontroldynamicsleftside.value);
        if ~isempty(implicitnonlinearcontrolindex.value)
            for ii=1:numterm
                if ii==1
                    if numterm>1
                        algebraic.term{ii}=['[' optimalcontroldynamicsleftside.value{ii} '; ...'];
                    else
                        algebraic.term{ii}=['[' optimalcontroldynamicsleftside.value{ii}];
                    end
                elseif ii<numterm
                    algebraic.term{ii}=[optimalcontroldynamicsleftside.value{ii} '; ...'];
                else
                    algebraic.term{ii}=optimalcontroldynamicsleftside.value{ii};
                end
            end
            algebraic.term{ii}=[algebraic.term{ii} ']'];
        else
            algebraic.term=[];
        end

    case 'optimalcontroldynamics'
        arcidentifier=field2arcidentifier(arcfield);
        optimalcontroldynamicsleftside=retrieveimpulsemodelinformation(ocStruct,'optimalcontroldynamics',arcidentifier,getsymkernel);
        %implicitnonlinearcontrolindex=retrieveimpulsemodelinformation(ocStruct,'implicitnonlinearcontrolindex',arcidentifier);
        implicitnonlinearcontrolindex=retrieveimpulsemodelinformation(ocStruct,'totalimplicitvariableindex',arcidentifier);
        numterm=numel(optimalcontroldynamicsleftside.value);
        if ~isempty(implicitnonlinearcontrolindex.value)
            for ii=1:numterm
                if ii==1
                    if numterm>1
                        algebraic.term{ii}=['[' optimalcontroldynamicsleftside.value{ii} '; ...'];
                    else
                        algebraic.term{ii}=optimalcontroldynamicsleftside.value{ii};
                    end
                elseif ii<numterm
                    algebraic.term{ii}=[optimalcontroldynamicsleftside.value{ii} '; ...'];
                else
                    algebraic.term{ii}=optimalcontroldynamicsleftside.value{ii};
                end
            end
            if numterm>1
                algebraic.term{ii}=[algebraic.term{ii} ']'];
            end
        else
            algebraic.term=[];
        end

    case 'symbolicoptimalcontroldynamics'
        arcidentifier=field2arcidentifier(arcfield);
        optimalcontroldynamicsleftside=retrieveimpulsemodelinformation(ocStruct,'optimalcontroldynamics',arcidentifier,getsymkernel);
        %implicitnonlinearcontrolindex=retrieveimpulsemodelinformation(ocStruct,'implicitnonlinearcontrolindex',arcidentifier);
        implicitnonlinearcontrolindex=retrieveimpulsemodelinformation(ocStruct,'totalimplicitvariableindex',arcidentifier);
        numterm=numel(optimalcontroldynamicsleftside.value);
        if ~isempty(implicitnonlinearcontrolindex.value)
            for ii=1:numterm
                if ii==1
                    if numterm>1
                        algebraic.term{ii}=['sym([''[' optimalcontroldynamicsleftside.value{ii} ';'' ...'];
                    else
                        algebraic.term{ii}=['sym(''' optimalcontroldynamicsleftside.value{ii} ''')'];
                    end
                elseif ii<numterm
                    algebraic.term{ii}=['''' optimalcontroldynamicsleftside.value{ii} ';'' ...'];
                else
                    algebraic.term{ii}=['''' optimalcontroldynamicsleftside.value{ii} ']''])'];
                end
            end
            if numterm>1
                algebraic.term{ii}=[algebraic.term{ii} ']'];
            end
        else
            algebraic.term=[];
        end

    case 'invoptimalcontroldynamicsleftside'
        arcidentifier=field2arcidentifier(arcfield);
        optimalcontroldynamicsleftside=retrieveimpulsemodelinformation(ocStruct,'invoptimalcontroldynamicsleftside',arcidentifier,getsymkernel);
        implicitnonlinearcontrolindex=retrieveimpulsemodelinformation(ocStruct,'implicitnonlinearcontrolindex',arcidentifier);
        numterm=numel(optimalcontroldynamicsleftside.value);
        if ~isempty(implicitnonlinearcontrolindex.value)
            for ii=1:numterm
                if ii==1
                    if numterm>1
                        algebraic.term{ii}=['[' optimalcontroldynamicsleftside.value{ii} '; ...'];
                    else
                        algebraic.term{ii}=['[' optimalcontroldynamicsleftside.value{ii}];
                    end
                elseif ii<numterm
                    algebraic.term{ii}=[optimalcontroldynamicsleftside.value{ii} '; ...'];
                else
                    algebraic.term{ii}=optimalcontroldynamicsleftside.value{ii};
                end
            end
            algebraic.term{ii}=[algebraic.term{ii} ']'];
        else
            algebraic.term=[];
        end

    case 'tensoroptimalcontroldynamicsleftside'
        arcidentifier=field2arcidentifier(arcfield);
        tensoroptimalcontroldynamicsleftside=retrieveimpulsemodelinformation(ocStruct,'tensoroptimalcontroldynamicsleftside',arcidentifier,getsymkernel);
        implicitnonlinearcontrolindex=retrieveimpulsemodelinformation(ocStruct,'implicitnonlinearcontrolindex',arcidentifier);
        nummat=numel(tensoroptimalcontroldynamicsleftside.value);
        if ~isempty(implicitnonlinearcontrolindex.value)
            counter=0;
            for jj=1:nummat
                numterm=numel(tensoroptimalcontroldynamicsleftside.value{jj});
                for ii=1:numterm
                    counter=counter+1;
                    if ii==1
                        if numterm>1
                            algebraic.term{counter}=['out(:,:,' num2str(jj) ')=[' tensoroptimalcontroldynamicsleftside.value{jj}{ii} '; ...'];
                        else
                            algebraic.term{counter}=['out(:,:,' num2str(jj) ')=[' tensoroptimalcontroldynamicsleftside.value{jj}{ii}];
                        end
                    elseif ii<numterm
                        algebraic.term{counter}=[tensoroptimalcontroldynamicsleftside.value{jj}{ii} '; ...'];
                    else
                        algebraic.term{counter}=tensoroptimalcontroldynamicsleftside.value{jj}{ii};
                    end
                end
                algebraic.term{counter}=[algebraic.term{counter} '];'];
            end
        else
            algebraic.term=[];
        end


    case 'jacobianoptimalcontroldynamicsrightside'
        arcidentifier=field2arcidentifier(arcfield);
        jacobianoptimalcontroldynamicsrightside=retrieveimpulsemodelinformation(ocStruct,'jacobianoptimalcontroldynamicsrightside',arcidentifier,getsymkernel);
        implicitnonlinearcontrolindex=retrieveimpulsemodelinformation(ocStruct,'implicitnonlinearcontrolindex',arcidentifier);
        numterm=numel(jacobianoptimalcontroldynamicsrightside.value);
        if ~isempty(implicitnonlinearcontrolindex.value)
            for ii=1:numterm
                if ii==1
                    if numterm>1
                        algebraic.term{ii}=['[' jacobianoptimalcontroldynamicsrightside.value{ii} '; ...'];
                    else
                        algebraic.term{ii}=['[' jacobianoptimalcontroldynamicsrightside.value{ii}];
                    end
                elseif ii<numterm
                    algebraic.term{ii}=[jacobianoptimalcontroldynamicsrightside.value{ii} '; ...'];
                else
                    algebraic.term{ii}=jacobianoptimalcontroldynamicsrightside.value{ii};
                end
            end
            algebraic.term{ii}=[algebraic.term{ii} ']'];
        else
            algebraic.term=[];
        end

    case 'pontryaginfunctionDuDX'
        pontryaginfunctionDuDX=retrieveimpulsemodelinformation(ocStruct,'pontryaginfunctionDuDX','','');
        for ii=1:numel(pontryaginfunctionDuDX.value)
            if ii==1
                if numel(pontryaginfunctionDuDX.value)>1
                    algebraic.term{ii}=['[' pontryaginfunctionDuDX.value{ii} '; ...'];
                else
                    algebraic.term{ii}=['[' pontryaginfunctionDuDX.value{ii}];
                end
            elseif ii<numel(pontryaginfunctionDuDX.value)
                algebraic.term{ii}=[pontryaginfunctionDuDX.value{ii} '; ...'];
            else
                algebraic.term{ii}=pontryaginfunctionDuDX.value{ii};
            end
        end
        algebraic.term{ii}=[algebraic.term{ii} ']'];

    case 'pontryaginfunctionDlmmcDX'
        arcidentifier=field2arcidentifier(arcfield);
        pontryaginfunctionDlmmcDX=retrieveimpulsemodelinformation(ocStruct,'pontryaginfunctionDlmmcDX',arcidentifier,getsymkernel);
        if ~isempty(pontryaginfunctionDlmmcDX.value)
            for ii=1:numel(pontryaginfunctionDlmmcDX.value)
                if ii==1
                    if numel(pontryaginfunctionDlmmcDX.value)>1
                        algebraic.term{ii}=['[' pontryaginfunctionDlmmcDX.value{ii} '; ...'];
                    else
                        algebraic.term{ii}=['[' pontryaginfunctionDlmmcDX.value{ii}];
                    end
                elseif ii<numel(pontryaginfunctionDlmmcDX.value)
                    algebraic.term{ii}=[pontryaginfunctionDlmmcDX.value{ii} '; ...'];
                else
                    algebraic.term{ii}=pontryaginfunctionDlmmcDX.value{ii};
                end
            end
            algebraic.term{ii}=[algebraic.term{ii} ']'];
        else
            algebraic.term=[];
        end


    case 'canonicalsystemparameterjacobian'
        arcidentifier=field2arcidentifier(arcfield);
        canonicalsystemparameterjacobian=retrieveimpulsemodelinformation(ocStruct,'canonicalsystemparameterjacobian',arcidentifier,getsymkernel);
        numterm=numel(canonicalsystemparameterjacobian.value);
        for ii=1:numterm
            if ii==1
                algebraic.term{ii}=['[' canonicalsystemparameterjacobian.value{ii} '; ...'];
            elseif ii<numterm
                algebraic.term{ii}=[canonicalsystemparameterjacobian.value{ii} '; ...'];
            else
                algebraic.term{ii}=canonicalsystemparameterjacobian.value{ii};
            end
        end
        algebraic.term{ii}=[algebraic.term{ii} ']'];


    case 'statecontrolcanonicalsystemparameterjacobian'
        arcidentifier=field2arcidentifier(arcfield);
        statecontrolcanonicalsystemparameterjacobian=retrieveimpulsemodelinformation(ocStruct,'statecontrolcanonicalsystemparameterjacobian',arcidentifier,getsymkernel);
        numterm=numel(statecontrolcanonicalsystemparameterjacobian.value);
        for ii=1:numterm
            if ii==1
                algebraic.term{ii}=['[' statecontrolcanonicalsystemparameterjacobian.value{ii} '; ...'];
            elseif ii<numterm
                algebraic.term{ii}=[statecontrolcanonicalsystemparameterjacobian.value{ii} '; ...'];
            else
                algebraic.term{ii}=statecontrolcanonicalsystemparameterjacobian.value{ii};
            end
        end
        algebraic.term{ii}=[algebraic.term{ii} ']'];

    case 'canonicalsystemderivativetime'
        arcidentifier=field2arcidentifier(arcfield);
        canonicalsystemderivativetime=retrieveimpulsemodelinformation(ocStruct,'canonicalsystemderivativetime',arcidentifier,getsymkernel);
        numterm=numel(canonicalsystemderivativetime.value);
        for ii=1:numterm
            if ii==1
                algebraic.term{ii}=['[' canonicalsystemderivativetime.value{ii} '; ...'];
            elseif ii<numterm
                algebraic.term{ii}=[canonicalsystemderivativetime.value{ii} '; ...'];
            else
                algebraic.term{ii}=canonicalsystemderivativetime.value{ii};
            end
        end
        algebraic.term{ii}=[algebraic.term{ii} ']'];

    case 'canonicalsystemhessian'
        arcidentifier=field2arcidentifier(arcfield);
        canonicalsystemhessian=retrieveimpulsemodelinformation(ocStruct,'canonicalsystemhessian',arcidentifier,getsymkernel);
        nummatrix=numel(canonicalsystemhessian.value);
        numlines=numel(canonicalsystemhessian.value{1});
        counter=0;
        for ii=1:nummatrix
            for jj=1:numlines
                counter=counter+1;
                if jj==1
                    algebraic.term{counter}=['[' canonicalsystemhessian.value{ii}{jj} '; ...'];
                elseif jj<numlines
                    algebraic.term{counter}=[canonicalsystemhessian.value{ii}{jj} '; ...'];
                elseif jj==numlines
                    algebraic.term{counter}=canonicalsystemhessian.value{ii}{jj};
                end
            end
            algebraic.term{counter}=[algebraic.term{counter} ']'];
        end
        algebraic.info=numlines; % return number of matrices for Hessian

    case 'canonicalsystemtotalhessian'
        arcidentifier=field2arcidentifier(arcfield);
        canonicalsystemparameterhessian=retrieveimpulsemodelinformation(ocStruct,'canonicalsystemtotalhessian',arcidentifier,getsymkernel);
        nummatrix=numel(canonicalsystemparameterhessian.value);
        numlines=numel(canonicalsystemparameterhessian.value{1});
        counter=0;
        for ii=1:nummatrix
            for jj=1:numlines
                counter=counter+1;
                if jj==1
                    algebraic.term{counter}=['[' canonicalsystemparameterhessian.value{ii}{jj} '; ...'];
                elseif jj<numlines
                    algebraic.term{counter}=[canonicalsystemparameterhessian.value{ii}{jj} '; ...'];
                elseif jj==numlines
                    algebraic.term{counter}=canonicalsystemparameterhessian.value{ii}{jj};
                end
            end
            algebraic.term{counter}=[algebraic.term{counter} ']'];
        end
        algebraic.info=numlines; % return number of matrices for Hessian

    case 'parametervariables'
        parametername=retrieveimpulsemodelinformation(ocStruct,'parametername');
        algebraic.term=parametername.value;

    case 'discobjectivefunction'
        discobjectivefunction=retrieveimpulsemodelinformation(ocStruct,'discobjectivefunction');
        algebraic.term{1}=discobjectivefunction.value;

    case 'discobjectivefunctionderivativetime'
        arcidentifier=field2arcidentifier(arcfield);
        discobjectivefunctionderivativetime=retrieveimpulsemodelinformation(ocStruct,'discobjectivefunctionderivativetime',arcidentifier,getsymkernel);
        algebraic.term=discobjectivefunctionderivativetime.value;

    case 'symbolicdiscobjectivefunction'
        discobjectivefunction=retrieveimpulsemodelinformation(ocStruct,'discobjectivefunction');
        algebraic.term{1}=['sym(''' discobjectivefunction.value ''')'];

    case 'discobjectivefunctionjacobian'
        arcidentifier=field2arcidentifier(arcfield);
        discobjectivefunctionjacobian=retrieveimpulsemodelinformation(ocStruct,'discobjectivefunctionDX',arcidentifier,getsymkernel);
        numterm=numel(discobjectivefunctionjacobian.value);
        for ii=1:numterm
            if ii==1
                algebraic.term{ii}=['[' discobjectivefunctionjacobian.value{ii} ', ...'];
            elseif ii<numterm
                algebraic.term{ii}=[discobjectivefunctionjacobian.value{ii} ', ...'];
            else
                algebraic.term{ii}=discobjectivefunctionjacobian.value{ii};
            end
        end
        algebraic.term{ii}=[algebraic.term{ii} ']'];

    case 'discobjectivefunctionparameterjacobian'
        arcidentifier=field2arcidentifier(arcfield);
        discobjectivefunctionjacobian=retrieveimpulsemodelinformation(ocStruct,'discobjectivefunctionDP',arcidentifier,getsymkernel);
        numterm=numel(discobjectivefunctionjacobian.value);
        for ii=1:numterm
            if ii==1
                algebraic.term{ii}=['[' discobjectivefunctionjacobian.value{ii} ', ...'];
            elseif ii<numterm
                algebraic.term{ii}=[discobjectivefunctionjacobian.value{ii} ', ...'];
            else
                algebraic.term{ii}=discobjectivefunctionjacobian.value{ii};
            end
        end
        algebraic.term{ii}=[algebraic.term{ii} ']'];

    case 'discimpulseobjectivefunction'
        discobjectivefunction=retrieveimpulsemodelinformation(ocStruct,'discimpulseobjectivefunction');
        algebraic.term{1}=discobjectivefunction.value;

    case 'symbolicdiscimpulseobjectivefunction'
        discimpulseobjectivefunction=retrieveimpulsemodelinformation(ocStruct,'discimpulseobjectivefunction');
        algebraic.term{1}=['sym(''' discimpulseobjectivefunction.value ''')'];

    case 'hamiltonianfunction'
        hamiltonianfunction=retrieveimpulsemodelinformation(ocStruct,'hamiltonianfunction');
        algebraic.term{1}=hamiltonianfunction.value;

    case 'impulsehamiltonianfunction'
        hamiltonianfunction=retrieveimpulsemodelinformation(ocStruct,'impulsehamiltonianfunction');
        algebraic.term{1}=hamiltonianfunction.value;

    case 'pontryaginfunction'
        pontryaginfunction=retrieveimpulsemodelinformation(ocStruct,'pontryaginfunction');
        algebraic.term{1}=pontryaginfunction.value;

    case 'salvagevalue'
        pontryaginfunction=retrieveimpulsemodelinformation(ocStruct,'salvagevalue');
        algebraic.term{1}=pontryaginfunction.value;

    case 'symbolicsalvagevalue'
        salvagevaluefunction=retrieveimpulsemodelinformation(ocStruct,'salvagevalue');
        algebraic.term{1}=['sym(''' salvagevaluefunction.value ''')'];

    case 'symbolicpontryaginfunction'
        pontryaginfunction=retrieveimpulsemodelinformation(ocStruct,'pontryaginfunction');
        algebraic.term{1}=['sym(''' pontryaginfunction.value ''')'];

    case 'symbolicimpulsepontryaginfunction'
        impulsepontryaginfunction=retrieveimpulsemodelinformation(ocStruct,'impulsehamiltonianfunction');
        algebraic.term{1}=['sym(''' impulsepontryaginfunction.value ''')'];

    case 'pontryaginfunctionDx'
        pontryaginfunctionDx=retrieveimpulsemodelinformation(ocStruct,'pontryaginfunctionDx');
        numterm=numel(pontryaginfunctionDx.value);
        for ii=1:numterm
            if ii==1
                algebraic.term{ii}=['[' pontryaginfunctionDx.value{ii} '; ...'];
            elseif ii<numterm
                algebraic.term{ii}=[pontryaginfunctionDx.value{ii} '; ...'];
            else
                algebraic.term{ii}=pontryaginfunctionDx.value{ii};
            end
        end
        algebraic.term{ii}=[algebraic.term{ii} ']'];

    case 'pontryaginfunctionDu2'
        pontryaginfunctionDu2=retrieveimpulsemodelinformation(ocStruct,'pontryaginfunctionDu2');
        numterm=numel(pontryaginfunctionDu2.value);
        for ii=1:numterm
            if ii==1
                algebraic.term{ii}=['[' pontryaginfunctionDu2.value{ii} '; ...'];
            elseif ii<numterm
                algebraic.term{ii}=[pontryaginfunctionDu2.value{ii} '; ...'];
            else
                algebraic.term{ii}=pontryaginfunctionDu2.value{ii};
            end
        end
        algebraic.term{ii}=[algebraic.term{ii} ']'];

    case 'pontryaginfunctionDX'
        arcidentifier=field2arcidentifier(arcfield);
        pontryaginfunctionDx=retrieveimpulsemodelinformation(ocStruct,'pontryaginfunctionDX',arcidentifier,getsymkernel);
        numterm=numel(pontryaginfunctionDx.value);
        for ii=1:numterm
            if ii==1
                algebraic.term{ii}=['[' pontryaginfunctionDx.value{ii} ', ...'];
            elseif ii<numterm
                algebraic.term{ii}=[pontryaginfunctionDx.value{ii} ', ...'];
            else
                algebraic.term{ii}=pontryaginfunctionDx.value{ii};
            end
        end
        algebraic.term{ii}=[algebraic.term{ii} ']'];

    case 'impulsepontryaginfunctionDtau'
        impulsepontryaginfunctionDtau=retrieveimpulsemodelinformation(ocStruct,'impulsepontryaginfunctionDtau');
        numterm=numel(impulsepontryaginfunctionDtau.value);
        for ii=1:numterm
            if ii==1 && numterm>1
                algebraic.term{ii}=['[' impulsepontryaginfunctionDtau.value{ii} '; ...'];
            elseif ii<numterm
                algebraic.term{ii}=[impulsepontryaginfunctionDtau.value{ii} '; ...'];
            else
                algebraic.term{ii}=impulsepontryaginfunctionDtau.value{ii};
            end
        end
        if numterm>1
            algebraic.term{ii}=[algebraic.term{ii} ']'];
        end

    case 'inequalitystateconstrainttimederivative'
        inequalitystateconstraint=retrieveimpulsemodelinformation(ocStruct,'inequalitystateconstraint');
        inequalitystateconstrainttimederivative=retrieveimpulsemodelinformation(ocStruct,'inequalitystateconstrainttimederivative');
        inequalitystateconstraintorder=retrieveimpulsemodelinformation(ocStruct,'inequalitystateconstraintorder');
        inequalitystateconstraintnum=retrieveimpulsemodelinformation(ocStruct,'inequalitystateconstraintnum');
        algebraic.term{1}='';
        for jj=1:inequalitystateconstraintnum.value
            algebraic.term{1}=[algebraic.term{1} inequalitystateconstraint.value{jj} ';'];
        end
        algebraic.term{1}(end)=[];
        for ii=(1:inequalitystateconstraintorder.value)+1
            algebraic.term{ii}='';
            for jj=1:inequalitystateconstraintnum.value
                algebraic.term{ii}=[algebraic.term{ii} inequalitystateconstrainttimederivative.value{jj}{ii-1} ';'];
            end
            algebraic.term{ii}(end)='';
        end


    case 'hamiltonianDu'
        hamiltonianDu=retrieveimpulsemodelinformation(ocStruct,'hamiltonianDu','',getsymkernel);
        numterm=numel(hamiltonianDu.value);
        for ii=1:numterm
            if ii==1
                if numterm>1
                    algebraic.term{ii}=['[' hamiltonianDu.value{ii} '; ...'];
                else
                    algebraic.term{ii}=[hamiltonianDu.value{ii}];
                end
            elseif ii<numterm
                algebraic.term{ii}=[hamiltonianDu.value{ii} '; ...'];
            else
                algebraic.term{ii}=hamiltonianDu.value{ii};
            end
        end
        if numterm>1
            algebraic.term{ii}=[algebraic.term{ii} ']'];
        end


    case 'DsalvagevalueDx'
        DsalvagevalueDx=retrieveimpulsemodelinformation(ocStruct,'DsalvagevalueDx','',getsymkernel);
        numterm=numel(DsalvagevalueDx.value);
        for ii=1:numterm
            if ii==1 && numterm>1
                algebraic.term{ii}=['[' DsalvagevalueDx.value{ii} '; ...'];
            elseif ii<numterm
                algebraic.term{ii}=[DsalvagevalueDx.value{ii} '; ...'];
            else
                algebraic.term{ii}=DsalvagevalueDx.value{ii};
            end
        end
        if numterm>1
            algebraic.term{ii}=[algebraic.term{ii} ']'];
        end
    case 'statedynamics'
        statedynamics=retrieveimpulsemodelinformation(ocStruct,'statedynamics','','');
        for ii=1:numel(statedynamics.value)
            if ii==1
                if numel(statedynamics.value)>1
                    algebraic.term{ii}=['[' statedynamics.value{ii} '; ...'];
                else
                    algebraic.term{ii}=[statedynamics.value{ii}];
                end
            elseif ii<numel(statedynamics.value)
                algebraic.term{ii}=[statedynamics.value{ii} '; ...'];
            else
                if numel(statedynamics.value)>1
                    algebraic.term{ii}=[statedynamics.value{ii} ']'];
                end
            end
        end

    case 'stateevent'
        stateevent=retrieveimpulsemodelinformation(ocStruct,'stateevent','','');
        for ii=1:numel(stateevent.value)
            if ii==1
                if numel(stateevent.value)>1
                    algebraic.term{ii}=['[' stateevent.value{ii} '; ...'];
                else
                    algebraic.term{ii}=[stateevent.value{ii}];
                end
            elseif ii<numel(stateevent.value)
                algebraic.term{ii}=[stateevent.value{ii} '; ...'];
            else
                if numel(stateevent.value)>1
                    algebraic.term{ii}=[stateevent.value{ii} ']'];
                end
            end
        end

    case 'adjointevent'
        adjointevent=retrieveimpulsemodelinformation(ocStruct,'adjointevent','','');
        for ii=1:numel(adjointevent.value)
            if ii==1
                if numel(adjointevent.value)>1
                    algebraic.term{ii}=['[' adjointevent.value{ii} '; ...'];
                else
                    algebraic.term{ii}=[adjointevent.value{ii}];
                end
            elseif ii<numel(adjointevent.value)
                algebraic.term{ii}=[adjointevent.value{ii} '; ...'];
            else
                if numel(adjointevent.value)>1
                    algebraic.term{ii}=[adjointevent.value{ii} ']'];
                end
            end
        end

    case 'adjointsystem'
        adjointsystem=retrieveimpulsemodelinformation(ocStruct,'adjointsystem','0','');
        for ii=1:numel(adjointsystem.value)
            if ii==1
                if numel(adjointsystem.value)>1
                    algebraic.term{ii}=['[' adjointsystem.value{ii} '; ...'];
                else
                    algebraic.term{ii}=[adjointsystem.value{ii}];
                end
            elseif ii<numel(adjointsystem.value)
                algebraic.term{ii}=[adjointsystem.value{ii} '; ...'];
            else
                if numel(adjointsystem.value)>1
                    algebraic.term{ii}=[adjointsystem.value{ii} ']'];
                end
            end
        end

    case 'objectivefunction'
        objectivefunction=retrieveimpulsemodelinformation(ocStruct,'objectiveintegrand');
        algebraic.term{1}=objectivefunction.value;

    case 'transversalitycondition'
        transversalitycondition=retrieveimpulsemodelinformation(ocStruct,'transversalitycondition','',getsymkernel);
        for ii=1:numel(transversalitycondition.value)
            if ii==1
                if numel(transversalitycondition.value)>1
                    algebraic.term{ii}=['[' transversalitycondition.value{ii} '; ...'];
                else
                    algebraic.term{ii}=[transversalitycondition.value{ii}];
                end
            elseif ii<numel(transversalitycondition.value)
                algebraic.term{ii}=[transversalitycondition.value{ii} '; ...'];
            else
                if numel(transversalitycondition.value)>1
                    algebraic.term{ii}=[transversalitycondition.value{ii} ']'];
                end
            end
        end
end

if vectflag
    if iscell(algebraic.term)
        for ii=1:numel(algebraic.term)
            algebraic.term{ii}=simplevectorize(algebraic.term{ii});
        end
    elseif ischar(algebraic.term)
        algebraic.term=simplevectorize(algebraic.term);
    end
end

function v = simplevectorize(v)
%VECTORIZE Vectorize expression.
%   VECTORIZE(S), when S is a string expression, inserts a '.' before
%   any '^', '*' or '/' in S.  The result is a character string.
%
%   VECTORIZE(FUN), when FUN is an INLINE function object, vectorizes the
%   formula for FUN.  The result is the vectorized version of the INLINE
%   function.
%
%   See also INLINE/FORMULA, INLINE.

%   Copyright 1984-2002 The MathWorks, Inc.
%   $Revision: 1.11 $  $Date: 2002/04/15 04:22:08 $

if isempty(v)
    v = [];
else
    for k = fliplr(find((v=='^') | (v=='*') | (v=='/')))
        v = [v(1:k-1) '.' v(k:end)];
    end
end