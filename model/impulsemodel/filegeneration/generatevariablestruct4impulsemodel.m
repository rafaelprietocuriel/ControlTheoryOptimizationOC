function varStruct=generatevariablestruct4impulsemodel(ocStruct,varStruct,varname,overwrite,varargin)
%
% a structure is generated which is then used to replace variables in
% template model files
% with overwrite 1 an already existing field of varStruct is overwritten
% the structure consists of a field given by the name of the variable. If
% the actual value of the variable depends on the arc the field string
% contains a cell array of strings and the field arcependent is set to one.
% multline 1 ... lines>1 are intended, e.g., CANONICALSYSTEMDYNAMICS
% multline 2 ...lines are not intended, e.g., PARAMETERVALUES
if nargin<=3
    overwrite=0;
end
if ~overwrite
    if isfield(varStruct,varname)
        return
    end
end
varStruct.(varname).arcdependent=0;
varStruct.(varname).arcidentifier=[];
varStruct.(varname).vectorize=0;
varStruct.(varname).type='';
varStruct.(varname).multline=0;

switch varname
    case 'INDEPENDENT'
        varStruct.(varname).string=getbasicname('independent');
        
    case 'JUMPTIME'
        varStruct.(varname).string=getbasicname('impulsetime');

    case 'ENDTIME'
        varStruct.(varname).string=getbasicname('endtime');

    case 'DEPENDENTVAR'
        varStruct.(varname).string=getbasicname('dependent');

    case 'IMPULSEDEPENDENTVAR'
        varStruct.(varname).string=getbasicname('impulsedependent');
        
    case 'INVIMPULSEDEPENDENTVAR'
        varStruct.(varname).string=['inv' getbasicname('impulsedependent')];
        
    case 'LDEPENDENTVAR'
         varStruct.(varname).string=[getbasicname('dependent') 'L'];

    case 'RDEPENDENTVAR'
         varStruct.(varname).string=[getbasicname('dependent') 'R'];

    case 'LRDEPENDENTVAR'
         varStruct.(varname).string=[getbasicname('dependent') 'LR'];

    case 'PARVAR'
        varStruct.(varname).string=getbasicname('parametervariables');

    case 'ARCVAR'
        varStruct.(varname).string=getbasicname('arcidentifiervar');

    case 'JUMPVAR'
        varStruct.(varname).string=getbasicname('jumpidentifiervar');

    case 'LRDEPENDENTVAR2L'
        depvar=getbasicname('dependent');
        lrdepvar=[getbasicname('dependent') 'LR'];

        varStruct.(varname).string=[depvar 'L=' lrdepvar '(:,1);'];
        varStruct.(varname).multline=0;
        
    case 'LRDEPENDENTVAR2R'
        depvar=getbasicname('dependent');
        lrdepvar=[getbasicname('dependent') 'LR'];

        varStruct.(varname).string=[depvar 'R=' lrdepvar '(:,2);'];
        varStruct.(varname).multline=0;

    case 'LRDEPENDENTVAR2IMPULSEDEPENDENTVAR'
        lrdepvar=[getbasicname('dependent') 'LR'];
        idepvar=getbasicname('impulsedependent');
        statenum=retrieveimpulsemodelinformation(ocStruct,'statenum');
        if statenum.value>1
            stateidx=['1:' num2str(statenum.value)];
            costateidx=[num2str(statenum.value+1) ':' num2str(2*statenum.value)];
        else
            stateidx='1';
            costateidx='2';
        end
        varStruct.(varname).string=[idepvar '=[' lrdepvar '(' stateidx ',1);' lrdepvar '(' costateidx ',2)];'];
        varStruct.(varname).multline=0;

    case 'LRDEPENDENTVAR2INVIMPULSEDEPENDENTVAR'
        lrdepvar=[getbasicname('dependent') 'LR'];
        invidepvar=['inv' getbasicname('impulsedependent')];
        statenum=retrieveimpulsemodelinformation(ocStruct,'statenum');
        if statenum.value>1
            stateidx=['1:' num2str(statenum.value)];
            costateidx=[num2str(statenum.value+1) ':' num2str(2*statenum.value)];
        else
            stateidx='1';
            costateidx='2';
        end
        varStruct.(varname).string=[invidepvar '=[' lrdepvar '(' stateidx ',2);' lrdepvar '(' costateidx ',1)];'];
        varStruct.(varname).multline=0;
  
    case 'IMPULSEDISCOUNTVARIABLE'
        sumdiscountratevariable=retrieveimpulsemodelinformation(ocStruct,'sumdiscountratevariable');
        varStruct.(varname).string=sumdiscountratevariable.value;
        varStruct.(varname).multline=0;
        
    case 'SUBSFLAG'
        varStruct.(varname).string='s';

    case 'LAGRANGEMULTCC'
        varStruct.(varname).string=getbasicname('lagrangemultcc');

    case 'LAGRANGEMULTSC'
        varStruct.(varname).string=getbasicname('lagrangemultsc');

    case 'CONTROL'
        varStruct.(varname).string=getbasicname('control');
        
    case 'IMPULSECONTROL'
        varStruct.(varname).string=getbasicname('impulsecontrol');

    case 'STATE'
        varStruct.(varname).string=getbasicname('state');

    case 'COSTATE'
        varStruct.(varname).string=getbasicname('costate');

    case 'EXOGENOUSFUNCTION'
        varStruct.(varname).string=getbasicname('exogenousfunction');

    case 'EXOGENOUSFUNCTIONDX'
        varStruct.(varname).string=['D1' getbasicname('exogenousfunction') '_Dx'];

    case 'GENERATIONDATE'
        varStruct.(varname).string=datestr(now);

    case 'MODELNAME'
        varStruct.(varname).string=modelname(ocStruct);

    case 'ARCIDENTIFIER'
        varStruct.(varname).arcdependent=1;
        arcidentifier=retrieveimpulsemodelinformation(ocStruct,'arcidentifier');
        varStruct.(varname).string=arcidentifier.value;

    case 'NUMBEROFARCS'
        arcnum=retrieveimpulsemodelinformation(ocStruct,'arcnum');
        varStruct.(varname).string=num2str(arcnum.value);

    case 'NUMBEROFJUMPIDS'
        arcnum=retrieveimpulsemodelinformation(ocStruct,'jumpidnum');
        varStruct.(varname).string=num2str(arcnum.value);

    case 'ARCARGUMENT'
        arcargument=retrieveimpulsemodelinformation(ocStruct,'argument');
        varStruct.(varname).string=mat2str(arcargument.value);

    case 'ARCCOLOR'
        arcnum=retrieveimpulsemodelinformation(ocStruct,'arcnum');
        colour=get(gca,'ColorOrder');
        if arcnum.value<=size(colour,1)
            varStruct.(varname).string=mat2str(colour(1:arcnum.value,:));
        else
            varStruct.(varname).string=mat2str([colour;rand(arcnum.value-size(colour,1),3)]);
        end

    case 'EDGES'
        arcargument=retrieveimpulsemodelinformation(ocStruct,'argument');
        arcnum=retrieveimpulsemodelinformation(ocStruct,'arcnum');
        if arcnum.value>1
            nc=nchoosek(arcargument.value,2);
            varStruct.(varname).string=mat2str([nc;nc(:,[2 1])]);
        else
            varStruct.(varname).string='[]';
        end
        
    case 'IMPULSECONTROLVALUESUBSTITUTION'
        icontrolnum=retrieveimpulsemodelinformation(ocStruct,'icontrolnum');
        icontrolname=retrieveimpulsemodelinformation(ocStruct,'icontrolname');
        indexstr='';
        if icontrolnum.value
            for ii=1:icontrolnum.value
                indexstr=[indexstr '''' icontrolname.value{ii} ''','];
            end
            coordstr=basename2vectorstring(getbasicname('impulsecontrol'),1:icontrolnum.value,'coord');
            indexstr(end)=[];
            varStruct.(varname).string=['out=subs(out,{' indexstr '},{' coordstr '});'];
        else
            varStruct.(varname).string='';
        end
    case 'CONTROLVALUESUBSTITUTION'
        controlnum=retrieveimpulsemodelinformation(ocStruct,'controlnum');
        controlname=retrieveimpulsemodelinformation(ocStruct,'controlname');
        inequalitycontrolconstraintnum=retrieveimpulsemodelinformation(ocStruct,'inequalitycontrolconstraintnum');
        lagrangemultipliercontrolname=retrieveimpulsemodelinformation(ocStruct,'lagrangemultipliercontrolname');
        indexstr='';
        if controlnum.value
            for ii=1:controlnum.value
                indexstr=[indexstr '''' controlname.value{ii} ''','];
            end
            coordstr=basename2vectorstring(getbasicname('control'),1:controlnum.value,'coord');
            if inequalitycontrolconstraintnum.value
                for ii=1:inequalitycontrolconstraintnum.value
                    indexstr=[indexstr '''' lagrangemultipliercontrolname.value{ii} ''','];
                end
                %indexstr=[indexstr ''',''' basename2vectorstring(getbasicname('lagrangemultcc'),1:inequalitycontrolconstraintnum.value,'index',''',''')];
                coordstr=[coordstr ',' basename2vectorstring(getbasicname('lagrangemultcc'),1:inequalitycontrolconstraintnum.value,'coord')];
            end
            indexstr(end)=[];
            varStruct.(varname).string=['out=subs(out,{' indexstr '},{' coordstr '});'];
        else
            varStruct.(varname).string='';
        end

    case 'COSTATEVALUESUBSTITUTION'
        costatenum=retrieveimpulsemodelinformation(ocStruct,'costatenum');
        costatename=retrieveimpulsemodelinformation(ocStruct,'costatename');
        indexstr='';
        for ii=1:costatenum.value
            indexstr=[indexstr '''' costatename.value{ii} ''','];
        end
        coordstr=basename2vectorstring(getbasicname('costate'),1:costatenum.value,'coord');
        indexstr(end)=[];
        varStruct.(varname).string=['out=subs(out,{' indexstr '},{' coordstr '});'];

    case 'STATEVALUESUBSTITUTION'
        statenum=retrieveimpulsemodelinformation(ocStruct,'statenum');
        statename=retrieveimpulsemodelinformation(ocStruct,'statename');
        indexstr='';
        for ii=1:statenum.value
            indexstr=[indexstr '''' statename.value{ii} ''','];
        end
        coordstr=basename2vectorstring(getbasicname('state'),1:statenum.value,'coord');
        indexstr(end)=[];
        varStruct.(varname).string=['out=subs(out,{' indexstr '},{' coordstr '});'];

    case 'CONSTRAINTNUM'
        inequalitycontrolconstraintnum=retrieveimpulsemodelinformation(ocStruct,'inequalitycontrolconstraintnum');
        inequalitystateconstraintnum=retrieveimpulsemodelinformation(ocStruct,'inequalitystateconstraintnum');
        varStruct.(varname).string=num2str(inequalitycontrolconstraintnum.value+inequalitystateconstraintnum.value);

    case 'CONTROLCONSTRAINTNUM'
        inequalitycontrolconstraintnum=retrieveimpulsemodelinformation(ocStruct,'inequalitycontrolconstraintnum');
        varStruct.(varname).string=num2str(inequalitycontrolconstraintnum.value);

    case 'ARCWITHSTATECONSTRAINT'
        arcwithinequalitystateconstraint=retrieveimpulsemodelinformation(ocStruct,'arcwithinequalitystateconstraint');
        varStruct.(varname).string=arcwithinequalitystateconstraint.value;

    case 'STATECONSTRAINTNUM'
        inequalitystateconstraintnum=retrieveimpulsemodelinformation(ocStruct,'inequalitystateconstraintnum');
        varStruct.(varname).string='';
        varStruct.(varname).string=num2str(inequalitystateconstraintnum.value);

    case 'STATECONSTRAINTORDER'
        inequalitystateconstraintorder=retrieveimpulsemodelinformation(ocStruct,'inequalitystateconstraintorder');
        try
            for ii=0:max(inequalitystateconstraintorder.value)
                varStruct.(varname).string{ii+1}=num2str(ii);
            end
            varStruct.(varname).arcdependent=1;
        catch
            %lasterr
            varStruct.(varname).string=[];
        end
       
    case 'STATECONSTRAINTBCR'
        %inequalitystateconstraintorder=retrieveimpulsemodelinformation(ocStruct,'inequalitystateconstraintorder');
        arcidentifier=retrieveimpulsemodelinformation(ocStruct,'arcidentifier');
        arcnum=retrieveimpulsemodelinformation(ocStruct,'arcnum');
        time=getbasicname('independent');
        dependent=getbasicname('dependent');
        parametervariables=getbasicname('parametervariables');
        for ii=1:arcnum.value
            fidx=retrieveimpulsemodelinformation(ocStruct,'nonzerolmscindex',arcidentifier.value{ii});
            if ~isempty(fidx.value)
                varStruct.(varname).string{ii}{1}=['DscDtN=' modelname(ocStruct) 'StateConstraintTotalTimeDerivative(' time ',' dependent 'r,' parametervariables ',edge(2),1);'];
                if length(fidx.value)==1
                    varStruct.(varname).string{ii}{2}=['out=DscDtN(' num2str(fidx.value) ');'];
                else
                    varStruct.(varname).string{ii}{2}=['out=DscDtN([' num2str(fidx.value) ']);'];
                end
            else
                varStruct.(varname).string{ii}{1}='out=[];';
            end
        end
        varStruct.(varname).arcdependent=1;
        varStruct.(varname).multline=1;

    case 'STATECONSTRAINTBCL'
        %inequalitystateconstraintorder=retrieveimpulsemodelinformation(ocStruct,'inequalitystateconstraintorder');
        arcidentifier=retrieveimpulsemodelinformation(ocStruct,'arcidentifier');
        arcnum=retrieveimpulsemodelinformation(ocStruct,'arcnum');
        time=getbasicname('independent');
        dependent=getbasicname('dependent');
        parametervariables=getbasicname('parametervariables');
        for ii=1:arcnum.value
            fidx=retrieveimpulsemodelinformation(ocStruct,'nonzerolmscindex',arcidentifier.value{ii});
            if ~isempty(fidx.value)
                varStruct.(varname).string{ii}{1}=['DscDtN=' modelname(ocStruct) 'StateConstraintTotalTimeDerivative(' time ',' dependent 'l,' parametervariables ',edge(1),0);'];
                if length(fidx.value)==1
                    varStruct.(varname).string{ii}{2}=['out=[out; DscDtN(' num2str(fidx.value) ')];'];
                else
                    varStruct.(varname).string{ii}{2}=['out=[out; DscDtN([' num2str(fidx.value) '])];'];
                end
                varStruct.(varname).string{ii}{3}=['DscDtN=' modelname(ocStruct) 'StateConstraintTotalTimeDerivative(' time ',' dependent 'l,' parametervariables ',edge(1),1);'];
                if length(fidx.value)==1
                    varStruct.(varname).string{ii}{4}=['out=[out; DscDtN(' num2str(fidx.value) ')];'];
                else
                    varStruct.(varname).string{ii}{4}=['out=[out; DscDtN([' num2str(fidx.value) '])];'];
                end
            end
        end
        varStruct.(varname).arcdependent=1;
        varStruct.(varname).multline=1;

    case 'STATECONSTRAINTTOTALTIMEDERIVATIVE'
        inequalitystateconstrainttimederivative=algebraicterm2string(ocStruct,'inequalitystateconstrainttimederivative',1,'');
        inequalitystateconstraintnum=retrieveimpulsemodelinformation(ocStruct,'inequalitystateconstraintnum');
        try
            for ii=1:numel(inequalitystateconstrainttimederivative.term)
                if inequalitystateconstraintnum.value==1
                    varStruct.(varname).string{ii}=['out=' inequalitystateconstrainttimederivative.term{ii}];
                else
                    varStruct.(varname).string{ii}=['out=[' inequalitystateconstrainttimederivative.term{ii} ']'];
                end
                varStruct.(varname).string{ii}=[varStruct.(varname).string{ii} ';'];
            end
            %end
            varStruct.(varname).multline=0;
            varStruct.(varname).vectorize=1;
            varStruct.(varname).arcdependent=1;
            varStruct.(varname).type='algterm';
        catch
            %lasterr
            varStruct.(varname).string=varname;
        end

    case 'ISAUTONOMOUS'
        autonomousflag=retrieveimpulsemodelinformation(ocStruct,'autonomous');
        varStruct.(varname).string=num2str(autonomousflag.value);

    case 'IMPLICITCONTROLS'
        implicitcontrolsflag=implicitcontrols(ocStruct);
        varStruct.(varname).string=num2str(implicitcontrolsflag);
 
    case 'NUMERICJACOBIAN'
        implicitcontrolsflag=implicitcontrols(ocStruct);
        varStruct.(varname).string=num2str(implicitcontrolsflag);
 
    case 'NUMERICHESSIAN'
        implicitcontrolsflag=implicitcontrols(ocStruct);
        varStruct.(varname).string=num2str(implicitcontrolsflag);

    case 'EXOGENOUSFUNCTIONNUM'
        exogenousfunctionnum=retrieveimpulsemodelinformation(ocStruct,'exogenousfunctionnum');
        varStruct.(varname).string=num2str(exogenousfunctionnum.value);

    case 'CONSTRAINT'
        try
            constraint=algebraicterm2string(ocStruct,'constraint',1);
            for jj=1:numel(constraint.term)
                if jj==1
                    varStruct.(varname).string{jj}=['out=[' constraint.term{jj}];
                    if numel(constraint.term)>1
                        varStruct.(varname).string{jj}=[varStruct.(varname).string{jj} '; ...'];
                    end
                elseif jj<numel(constraint.term)
                    varStruct.(varname).string{jj}=[constraint.term{jj} '; ...'];
                else
                    varStruct.(varname).string{jj}=constraint.term{jj};
                end
            end
            varStruct.(varname).string{jj}=[varStruct.(varname).string{jj} '];'];
            varStruct.(varname).multline=1;
            varStruct.(varname).vectorize=1;
            varStruct.(varname).arcdependent=0;
            varStruct.(varname).type='algterm';
        catch
            %lasterr
            varStruct.(varname).string=varname;
        end

    case 'STATECOSTATECOORD'
        statenum=retrieveimpulsemodelinformation(ocStruct,'statenum');
        costatenum=retrieveimpulsemodelinformation(ocStruct,'costatenum');
        varStruct.(varname).string=['1:' int2str(statenum.value+costatenum.value)];

    case 'CANONICALSYSTEMEQUATIONCOORD'
        canonicalsystemequationnum=retrieveimpulsemodelinformation(ocStruct,'canonicalsystemequationnum');
        varStruct.(varname).string=['1:' int2str(canonicalsystemequationnum.value)];

    case 'STATECOSTATENUM'
        statenum=retrieveimpulsemodelinformation(ocStruct,'statenum');
        costatenum=retrieveimpulsemodelinformation(ocStruct,'costatenum');
        varStruct.(varname).string=int2str(statenum.value+costatenum.value);

    case 'STATECOORD'
        statenum=retrieveimpulsemodelinformation(ocStruct,'statenum');
        varStruct.(varname).string=['1:' int2str(statenum.value)];

    case 'STATENUM'
        statenum=retrieveimpulsemodelinformation(ocStruct,'statenum');
        varStruct.(varname).string=int2str(statenum.value);

    case 'DEFAULTPLOTCOMMAND'
        depvarname=getbasicname('dependent');
        varStruct.(varname).string=['h=plot(' depvarname '(1,:),' depvarname '(2,:));'];

    case 'DEFAULTPLOTTHREECOMMAND'
        depvarname=getbasicname('dependent');
        varStruct.(varname).string=['h=plot3(' depvarname '(1,:),' depvarname '(2,:),' depvarname '(3,:));'];

    case 'DEFAULTPLOTINDIFFERENCECOMMAND'
        depvarname=getbasicname('dependent');
        varStruct.(varname).string{1}=['h=plot(' depvarname '(1,leftarcindex(1):rightarcindex(1)),' depvarname '(2,leftarcindex(1):rightarcindex(1)), ...'];
        varStruct.(varname).string{2}=['    ' depvarname '(1,leftarcindex(2):rightarcindex(2)),' depvarname '(2,leftarcindex(2):rightarcindex(2)));'];
        varStruct.(varname).multline=1;

    case 'DEFAULTPLOTTHREEINDIFFERENCECOMMAND'
        depvarname=getbasicname('dependent');
        varStruct.(varname).string{1}=['h=plot3(' depvarname '(1,leftarcindex(1):rightarcindex(1)),' depvarname '(2,leftarcindex(1):rightarcindex(1)),' depvarname '(3,leftarcindex(1):rightarcindex(1)), ...'];
        varStruct.(varname).string{2}=['    ' depvarname '(1,leftarcindex(2):rightarcindex(2)),' depvarname '(2,leftarcindex(2):rightarcindex(2)),' depvarname '(3,leftarcindex(2):rightarcindex(2)));'];
        varStruct.(varname).multline=1;

    case 'LATEXUSERDEPENDENTNAMEFIRST'
        statename=retrieveimpulsemodelinformation(ocStruct,'statename');
        [val1 val2 val3]=regexp(statename.value{1},'([0-9]+)\>','split');
        if ~isempty(val2)
            varStruct.(varname).string=['$' latex(sym(val1{1})) '_' statename.value{1}(val2:val3) '$'];
        else
            varStruct.(varname).string=['$' latex(sym(val1{1})) '$'];
        end
        %varStruct.(varname).string=['$' regexprep(statename.value{1},'([0-9]+)\>','_$1') '$'];

    case 'LATEXUSERDEPENDENTNAMESECOND'
        statenum=retrieveimpulsemodelinformation(ocStruct,'statenum');
        if statenum.value>1
            statename=retrieveimpulsemodelinformation(ocStruct,'statename');
            labelname=statename.value{2};
        else
            costatename=retrieveimpulsemodelinformation(ocStruct,'costatename');
            labelname=costatename.value{1};
        end
        [val1 val2 val3]=regexp(labelname,'([0-9]+)\>','split');
        if ~isempty(val2)
            varStruct.(varname).string=['$' latex(sym(val1{1})) '_' labelname(val2:val3) '$'];
        else
            varStruct.(varname).string=['$' latex(sym(val1{1})) '$'];
        end

    case 'LATEXUSERDEPENDENTNAMETHIRD'
        statenum=retrieveimpulsemodelinformation(ocStruct,'statenum');
        if statenum.value<=2
            labelname='';
        else
            statename=retrieveimpulsemodelinformation(ocStruct,'statename');
            labelname=statename.value{3};
            %varStruct.(varname).string=['$' regexprep(statename.value{2},'([0-9]+)\>','_$1') '$'];
        end

        [val1 val2 val3]=regexp(labelname,'([0-9]+)\>','split');
        if ~isempty(val2)
            varStruct.(varname).string=['$' latex(sym(val1{1})) '_' labelname(val2:val3) '$'];
        else
            varStruct.(varname).string=['$' latex(sym(val1{1})) '$'];
        end

    case 'EQUATIONNUM'
        try
            arcfield=getarcclass(ocStruct,'algebraicequation');
            for ii=1:numel(arcfield)
                arcidentifier=field2arcidentifier(arcfield{ii});
                equationnum=retrieveimpulsemodelinformation(ocStruct,'canonicalsystemequationnum',arcidentifier{1});
                varStruct.(varname).string{ii}=int2str(equationnum.value);
            end
            varStruct.(varname).arcdependent=1;
        catch
            %lasterr
            varStruct.(varname).string=varname;
        end
        
    case 'EQUATIONNUMIMPLICIT'
        try
            statenum=retrieveimpulsemodelinformation(ocStruct,'statenum');
            costatenum=retrieveimpulsemodelinformation(ocStruct,'costatenum');
            totalalgebraicequationnum=retrieveimpulsemodelinformation(ocStruct,'totalalgebraicequationnum');
            varStruct.(varname).string=int2str(statenum.value+costatenum.value+totalalgebraicequationnum.value);
            varStruct.(varname).arcdependent=0;
        catch
            %lasterr
            varStruct.(varname).string=varname;
        end

    case 'ALGEBRAICEQUATIONNUM'
        % returns
        try
            arcfield=getarcclass(ocStruct,'algebraicequation');
            algebraicequationnum=0;
            counter=1;
            while counter<=numel(arcfield)
                arcidentifier=field2arcidentifier(arcfield{counter});
                equationnum=retrieveimpulsemodelinformation(ocStruct,'algebraicequationnum',arcidentifier{1});
                if equationnum.value
                    algebraicequationnum=1;
                    break
                end
                counter=counter+1;
            end
            varStruct.(varname).string=int2str(algebraicequationnum);
            varStruct.(varname).arcdependent=0;
        catch
            %lasterr
            varStruct.(varname).string=varname;
        end

    case 'IMPLICITVARIABLENUM'
        % returns
        try
            arcfield=getarcclass(ocStruct,'algebraicequation');
            arcnum=retrieveimpulsemodelinformation(ocStruct,'arcnum');
            for ii=1:numel(arcfield)
                arcidentifier=field2arcidentifier(arcfield{ii});
                equationvariablenameimplicit=retrieveimpulsemodelinformation(ocStruct,'equationvariablenameimplicit',arcidentifier);
                varStruct.(varname).string{ii}=int2str(length(equationvariablenameimplicit.value));
            end
            if arcnum.value>1
                varStruct.(varname).arcdependent=1;
            else
                varStruct.(varname).arcdependent=0;
                varStruct.(varname).string=varStruct.(varname).string{1};
            end
            varStruct.(varname).multline=0;
        catch
            %lasterr
            varStruct.(varname).string=varname;
        end

    case 'TOTALIMPLICITVARIABLENUM'
        % returns
        try
            arcfield=getarcclass(ocStruct,'algebraicequation');
            arcnum=retrieveimpulsemodelinformation(ocStruct,'arcnum');
            for ii=1:numel(arcfield)
                totalalgebraicequationnum=retrieveimpulsemodelinformation(ocStruct,'totalalgebraicequationnum');
                varStruct.(varname).string{ii}=int2str(totalalgebraicequationnum.value);
            end
            if arcnum.value>1
                varStruct.(varname).arcdependent=1;
            else
                varStruct.(varname).arcdependent=0;
                varStruct.(varname).string=varStruct.(varname).string{1};
            end
            varStruct.(varname).multline=0;
        catch
            %lasterr
            varStruct.(varname).string=varname;
        end

    case 'IMPLICITNUM'
        % returns
        try
            implicitnonlinearcontrolindex=retrieveimpulsemodelinformation(ocStruct,'totalimplicitvariableindex');
            varStruct.(varname).string=int2str(length(implicitnonlinearcontrolindex.value));
            varStruct.(varname).arcdependent=0;
            varStruct.(varname).multline=0;
        catch
            %lasterr
            varStruct.(varname).string=varname;
        end

    case 'ARCALGEBRAICEQUATIONNUM'
        % returns
        try
            arcfield=getarcclass(ocStruct,'algebraicequation');
            arcnum=retrieveimpulsemodelinformation(ocStruct,'arcnum');
            for ii=1:numel(arcfield)
                arcidentifier=field2arcidentifier(arcfield{ii});
                equationnum=retrieveimpulsemodelinformation(ocStruct,'algebraicequationnum',arcidentifier);
                varStruct.(varname).string{ii}=int2str(equationnum.value);
            end
            if arcnum.value>1
                varStruct.(varname).arcdependent=1;
            else
                varStruct.(varname).arcdependent=0;
                varStruct.(varname).string=varStruct.(varname).string{1};
            end
            varStruct.(varname).multline=0;
        catch
            %lasterr
            varStruct.(varname).string=varname;
        end

    case 'ARCALGEBRAICEQUATIONNUMIMPLICIT'
        % returns
        try
            arcfield=getarcclass(ocStruct,'algebraicequation');
            arcnum=retrieveimpulsemodelinformation(ocStruct,'arcnum');
            for ii=1:numel(arcfield)
                arcidentifier=field2arcidentifier(arcfield{ii});
                equationnum=retrieveimpulsemodelinformation(ocStruct,'totalalgebraicequationnum',arcidentifier);
                varStruct.(varname).string{ii}=int2str(equationnum.value);
            end
            if arcnum.value>1
                varStruct.(varname).arcdependent=1;
            else
                varStruct.(varname).arcdependent=0;
                varStruct.(varname).string=varStruct.(varname).string{1};
            end
            varStruct.(varname).multline=0;
        catch
            %lasterr
            varStruct.(varname).string=varname;
        end

    case 'COSTATECOORD'
        statenum=retrieveimpulsemodelinformation(ocStruct,'statenum');
        costatenum=retrieveimpulsemodelinformation(ocStruct,'costatenum');
        varStruct.(varname).string=[int2str(statenum.value+1) ':' int2str(statenum.value+costatenum.value)];

    case 'IMPLICITCONTROLCOORD'
        varStruct.(varname).arcdependent=1;
        arcidentifier=retrieveimpulsemodelinformation(ocStruct,'arcidentifier');
        arcnum=retrieveimpulsemodelinformation(ocStruct,'arcnum');
        statenum=retrieveimpulsemodelinformation(ocStruct,'statenum');
        costatenum=retrieveimpulsemodelinformation(ocStruct,'costatenum');
        for ii=1:arcnum.value
            fidx=retrieveimpulsemodelinformation(ocStruct,'implicitnonlinearcontrolindex',arcidentifier.value{ii});
            if isempty(fidx.value)
                varStruct.(varname).string{ii}='[]';
            else
                if numel(fidx.value)>1
                    varStruct.(varname).string{ii}=['[' int2str(statenum.value+costatenum.value+(1:numel(fidx.value))) ']'];
                else
                    varStruct.(varname).string{ii}=int2str(statenum.value+costatenum.value+1);
                end
            end
        end

    case 'IMPLICITCONTROLINDEX'
        varStruct.(varname).arcdependent=1;
        arcidentifier=retrieveimpulsemodelinformation(ocStruct,'arcidentifier');
        arcnum=retrieveimpulsemodelinformation(ocStruct,'arcnum');
        for ii=1:arcnum.value
            fidx=retrieveimpulsemodelinformation(ocStruct,'implicitnonlinearcontrolindex',arcidentifier.value{ii});
            if isempty(fidx.value)
                varStruct.(varname).string{ii}='[]';
            else
                if length(fidx.value)==1
                    varStruct.(varname).string{ii}=int2str(fidx.value);
                else
                    varStruct.(varname).string{ii}=['[' int2str(fidx.value) ']'];
                end
            end
        end

    case 'OPTIMALCONTROLDYNAMICSINDEX'
        arcidentifier=retrieveimpulsemodelinformation(ocStruct,'arcidentifier');
        statenum=retrieveimpulsemodelinformation(ocStruct,'statenum');
        costatenum=retrieveimpulsemodelinformation(ocStruct,'costatenum');
        %totalimplicitnonlinearcontrolindex=retrieveimpulsemodelinformation(ocStruct,'totalimplicitvariableindex');
        arcnum=retrieveimpulsemodelinformation(ocStruct,'arcnum');
        for ii=1:arcnum.value
            fidx=retrieveimpulsemodelinformation(ocStruct,'implicitnonlinearcontrolindex',arcidentifier.value{ii});
            idx=zeros(1,numel(fidx.value));
            for jj=1:numel(fidx.value)
                idx(jj)=jj;%find(totalimplicitnonlinearcontrolindex.value==fidx.value(jj));
            end
            if ~isempty(fidx.value)
                if length(fidx.value)==1
                    varStruct.(varname).string{ii}=['DXDt(' int2str(statenum.value+costatenum.value+idx) ',:)=tmp(' int2str(fidx.value) ',:);'];
                else
                    varStruct.(varname).string{ii}=['DXDt([' int2str(statenum.value+costatenum.value+idx) '],:)=tmp([' int2str(fidx.value) '],:);'];
                end
            else
                varStruct.(varname).string{ii}='tmp=DXDt;';
            end
        end
        if arcnum.value>1
            varStruct.(varname).arcdependent=1;
        else
            varStruct.(varname).arcdependent=0;
            varStruct.(varname).string=varStruct.(varname).string{1};
        end
        varStruct.(varname).multline=0;
        varStruct.(varname).vectorize=0;

    case 'OPTIMALCONTROLDYNAMICSJACOBIANINDEX'
        varStruct.(varname).arcdependent=1;
        arcidentifier=retrieveimpulsemodelinformation(ocStruct,'arcidentifier');
        statenum=retrieveimpulsemodelinformation(ocStruct,'statenum');
        costatenum=retrieveimpulsemodelinformation(ocStruct,'costatenum');
        fidx=retrieveimpulsemodelinformation(ocStruct,'totalimplicitvariableindex');
        arcnum=retrieveimpulsemodelinformation(ocStruct,'arcnum');
        for ii=1:arcnum.value
            %fidx=retrieveimpulsemodelinformation(ocStruct,'implicitnonlinearcontrolindex',arcidentifier.value{ii});
            idx=zeros(1,numel(fidx.value));
            for jj=1:numel(fidx.value)
                idx(jj)=jj;%find(totalimplicitnonlinearcontrolindex.value==fidx.value(jj));
            end
            if ~isempty(fidx.value)
                if length(fidx.value)==1
                    varStruct.(varname).string{ii}=['DDXDtDX(' int2str(statenum.value+costatenum.value+idx) ',:)=tmp(' int2str(fidx.value) ',:);'];
                else
                    varStruct.(varname).string{ii}=['DDXDtDX([' int2str(statenum.value+costatenum.value+idx) '],:)=tmp([' int2str(fidx.value) '],:);'];
                end
            else
                varStruct.(varname).string{ii}='tmp=DDXDtDX;';
            end
        end
        if arcnum.value>1
            varStruct.(varname).arcdependent=1;
        else
            varStruct.(varname).arcdependent=0;
            varStruct.(varname).string=varStruct.(varname).string{1};
        end
        varStruct.(varname).multline=0;
        varStruct.(varname).vectorize=0;

    case 'ODEDIM'
        varStruct.(varname).arcdependent=1;
        odedim=retrieveimpulsemodelinformation(ocStruct,'odedim');
        arcnum=retrieveimpulsemodelinformation(ocStruct,'arcnum');
        for ii=1:arcnum.value
            varStruct.(varname).string{ii}=int2str(odedim.value);
        end

    case 'AEDIM'
        varStruct.(varname).arcdependent=1;
        arcnum=retrieveimpulsemodelinformation(ocStruct,'arcnum');
        arcidentifier=retrieveimpulsemodelinformation(ocStruct,'arcidentifier');
        for ii=1:arcnum.value
            algebraicequationnum=retrieveimpulsemodelinformation(ocStruct,'algebraicequationnum',arcidentifier.value{ii});
            varStruct.(varname).string{ii}=int2str(algebraicequationnum.value);
        end

    case 'DAEORDER'
        if ~isfield(varStruct,'AEDIM')
            varStruct=generatevariablestruct4standardmodel(ocStruct,varStruct,'AEDIM');
        end
        varStruct.(varname).arcdependent=1;
        arcnum=retrieveimpulsemodelinformation(ocStruct,'arcnum');
        statenum=retrieveimpulsemodelinformation(ocStruct,'statenum');
        costatenum=retrieveimpulsemodelinformation(ocStruct,'costatenum');
        for ii=1:arcnum.value
            order=['[' int2str(ones(1,statenum.value+costatenum.value))];
            aedim=str2double(varStruct.AEDIM.string{ii});
            if aedim>0
                order=[order '  ' int2str(zeros(1, aedim))];
            end
            order=[order ']'];
            varStruct.(varname).string{ii}=order;
        end
        
    case 'PARAMETERVALUES'
        % the strings used for assigning parametervalues from the par
        % argument in functions to the spcific paraemter variable names,
        % e.g. r=pararg(1); ...
        varStruct.(varname).multline=2;
        basename=getbasicname('parametervariables');
        parvar=algebraicterm2string(ocStruct,'parametervariables');
        parvarexpr=retrieveimpulsemodelinformation(ocStruct,'parametervalueexpression');
        for ii=1:numel(parvar.term)
            if isempty(parvarexpr.value{ii})
            varStruct.(varname).string{ii}=[parvar.term{ii} '=' basename '(' num2str(ii) ');'];
            else
                expr=parvarexpr.value{ii};
                for jj=1:numel(parvar.term)
                    if findstr(expr,parvar.term{jj})
                        expr=regexprep(expr,['\<' parvar.term{jj} '\>'],[basename '(' num2str(jj) ')']);
                    end
                end
                varStruct.(varname).string{ii}=[parvar.term{ii} '=' expr ';'];
            end
        end

    case 'ZEROSNTIMESNUMEQ'
        % returns a zero matrix where number of rows = number of states and
        % number of columns = number of equations of the canonical system
        try
            arcfield=getarcclass(ocStruct,'algebraicequation');
            statenum=retrieveimpulsemodelinformation(ocStruct,'statenum');
            for ii=1:numel(arcfield)
                arcidentifier=field2arcidentifier(arcfield{ii});
                equationnum=retrieveimpulsemodelinformation(ocStruct,'canonicalsystemequationnum',arcidentifier{1});
                zerovec=sprintf('%d,',zeros(1,equationnum.value));
                zerovec(end)=[];
                for jj=1:statenum.value
                    if jj==1
                        if statenum.value>1
                            varStruct.(varname).string{ii}{jj}=['[' zerovec '; ...'];
                        else
                            varStruct.(varname).string{ii}{jj}=['[' zerovec];
                        end
                    elseif jj<statenum.value
                        varStruct.(varname).string{ii}{jj}=[zerovec '; ...'];
                    else
                        varStruct.(varname).string{ii}{jj}=zerovec;
                    end
                end
                varStruct.(varname).string{ii}{jj}=[varStruct.(varname).string{ii}{jj} '];'];
            end
            varStruct.(varname).multline=2;
            varStruct.(varname).arcdependent=1;
        catch
            %lasterr
            varStruct.(varname).string=varname;
        end

    case 'ZEROSNUMEQ'
        % returns a zero matrix where number of rows = number of states and
        % number of columns = number of equations of the canonical system
        try
            arcfield=getarcclass(ocStruct,'algebraicequation');
            for ii=1:numel(arcfield)
                arcidentifier=field2arcidentifier(arcfield{ii});
                equationnum=retrieveimpulsemodelinformation(ocStruct,'canonicalsystemequationnum',arcidentifier{1});
                zerovec=sprintf('%d,',zeros(1,equationnum.value));
                zerovec(end)=[];
                varStruct.(varname).string{ii}=['[' zerovec ']'];
            end
            varStruct.(varname).multline=0;
            varStruct.(varname).arcdependent=1;
        catch
            %lasterr
            varStruct.(varname).string=varname;
        end

    case 'CANONICALSYSTEM'
        try
            canonicalsystem=algebraicterm2string(ocStruct,'canonicalsystem',1,'');
            for jj=1:numel(canonicalsystem.term)
                if jj==1
                    varStruct.(varname).string{jj}=canonicalsystem.term{jj};
                else
                    varStruct.(varname).string{jj}=canonicalsystem.term{jj};
                end
            end
            varStruct.(varname).string{jj}=[varStruct.(varname).string{jj} ';'];
            %end
            varStruct.(varname).multline=1;
            varStruct.(varname).vectorize=1;
            varStruct.(varname).arcdependent=0;
            varStruct.(varname).type='algterm';
        catch
            %lasterr
            varStruct.(varname).string=varname;
        end

    case 'EXOGENOUSFUNCTIONTERM'
        try
            exogenousfunctionterm=algebraicterm2string(ocStruct,'exogenousfunctionterm',1,'');
            for jj=1:numel(exogenousfunctionterm.term)
                if jj==1
                    varStruct.(varname).string{jj}=exogenousfunctionterm.term{jj};
                else
                    varStruct.(varname).string{jj}=exogenousfunctionterm.term{jj};
                end
            end
            varStruct.(varname).string{jj}=[varStruct.(varname).string{jj} ';'];
            %end
            varStruct.(varname).multline=1;
            varStruct.(varname).vectorize=1;
            varStruct.(varname).arcdependent=0;
            varStruct.(varname).type='algterm';
        catch
            %lasterr
            varStruct.(varname).string=varname;
        end

    case 'CANONICALSYSTEMDYNAMICS'
        try
            canonicalsystem=algebraicterm2string(ocStruct,'canonicalsystem',1,'');
            for jj=1:numel(canonicalsystem.term)
                if jj==1
                    varStruct.(varname).string{jj}=['out=' canonicalsystem.term{jj}];
                else
                    varStruct.(varname).string{jj}=canonicalsystem.term{jj};
                end
            end
            varStruct.(varname).string{jj}=[varStruct.(varname).string{jj} ';'];
            %end
            varStruct.(varname).multline=1;
            varStruct.(varname).vectorize=1;
            varStruct.(varname).arcdependent=0;
            varStruct.(varname).type='algterm';
        catch
            %lasterr
            varStruct.(varname).string=varname;
        end

    case 'CANONICALSYSTEMALGEBRAICEQUATION'
        try
            arcfield=getarcclass(ocStruct,'algebraicequation');
            arcnum=retrieveimpulsemodelinformation(ocStruct,'arcnum');
            varStruct.(varname).arcidentifier=field2arcidentifier(arcfield);
            for ii=1:numel(arcfield)
                algebraicequation=algebraicterm2string(ocStruct,'algebraicequation',1,arcfield{ii});
                if ~isempty(algebraicequation.term)
                    for jj=1:numel(algebraicequation.term)
                        if jj==1
                            varStruct.(varname).string{ii}{jj}=['outae=' algebraicequation.term{jj}];
                        else
                            varStruct.(varname).string{ii}{jj}=algebraicequation.term{jj};
                        end
                    end
                    varStruct.(varname).string{ii}{jj}=[varStruct.(varname).string{ii}{jj} ';'];
                else
                    varStruct.(varname).string{ii}{1}='outae=[];';
                end
            end
            varStruct.(varname).multline=1;
            varStruct.(varname).vectorize=1;
            if arcnum.value>1
                varStruct.(varname).arcdependent=1;
            else
                varStruct.(varname).arcdependent=0;
                varStruct.(varname).string=varStruct.(varname).string{1};
            end
            varStruct.(varname).type='algterm';
        catch
            %lasterr
            varStruct.(varname).string=varname;
        end

    case 'ALGEBRAICEQUATIONIMPLICIT'
        try
            arcfield=getarcclass(ocStruct,'algebraicequation');
            arcnum=retrieveimpulsemodelinformation(ocStruct,'arcnum');
            varStruct.(varname).arcidentifier=field2arcidentifier(arcfield);
            for ii=1:numel(arcfield)
                algebraicequationimplicit=algebraicterm2string(ocStruct,'algebraicequationimplicit',1,arcfield{ii});
                if ~isempty(algebraicequationimplicit.term)
                    for jj=1:numel(algebraicequationimplicit.term)
                        if jj==1
                            varStruct.(varname).string{ii}{jj}=['outae=' algebraicequationimplicit.term{jj}];
                        else
                            varStruct.(varname).string{ii}{jj}=algebraicequationimplicit.term{jj};
                        end
                    end
                    varStruct.(varname).string{ii}{jj}=[varStruct.(varname).string{ii}{jj} ';'];
                else
                    varStruct.(varname).string{ii}{1}='outae=[];';
                end
            end
            varStruct.(varname).multline=1;
            varStruct.(varname).vectorize=1;
            if arcnum.value>1
                varStruct.(varname).arcdependent=3;
            else
                varStruct.(varname).arcdependent=0;
                varStruct.(varname).string=varStruct.(varname).string{1};
            end
            varStruct.(varname).type='algterm';
        catch
            %lasterr
            varStruct.(varname).string=varname;
        end


    case 'INVOPTIMALCONTROLDYNAMICSLEFTSIDE'
        try
            arcfield=getarcclass(ocStruct,'canonicalsystemjacobian');
            arcnum=retrieveimpulsemodelinformation(ocStruct,'arcnum');
            varStruct.(varname).arcidentifier=field2arcidentifier(arcfield);
            for ii=1:numel(arcfield)
                optimalcontroldynamicsleftside=algebraicterm2string(ocStruct,'invoptimalcontroldynamicsleftside',0,arcfield{ii});
                if ~isempty(optimalcontroldynamicsleftside.term)
                    for jj=1:numel(optimalcontroldynamicsleftside.term)
                        if jj==1
                            varStruct.(varname).string{ii}{jj}=['out=' optimalcontroldynamicsleftside.term{jj}];
                        else
                            varStruct.(varname).string{ii}{jj}=optimalcontroldynamicsleftside.term{jj};
                        end
                    end
                    varStruct.(varname).string{ii}{jj}=[varStruct.(varname).string{ii}{jj} ';'];
                else
                    varStruct.(varname).string{ii}{1}='out=[];';
                end
            end
            varStruct.(varname).multline=1;
            if arcnum.value>1
                varStruct.(varname).arcdependent=1;
            else
                varStruct.(varname).arcdependent=0;
                varStruct.(varname).string=varStruct.(varname).string{1};
            end
            varStruct.(varname).type='algterm';
        catch
            %lasterr
            varStruct.(varname).string=varname;
        end

    case 'OPTIMALCONTROLDYNAMICSLEFTSIDE'
        try
            arcfield=getarcclass(ocStruct,'canonicalsystemjacobian');
            arcnum=retrieveimpulsemodelinformation(ocStruct,'arcnum');
            varStruct.(varname).arcidentifier=field2arcidentifier(arcfield);
            for ii=1:numel(arcfield)
                optimalcontroldynamicsleftside=algebraicterm2string(ocStruct,'optimalcontroldynamicsleftside',0,arcfield{ii});
                if ~isempty(optimalcontroldynamicsleftside.term)
                    for jj=1:numel(optimalcontroldynamicsleftside.term)
                        if jj==1
                            varStruct.(varname).string{ii}{jj}=['out=' optimalcontroldynamicsleftside.term{jj}];
                        else
                            varStruct.(varname).string{ii}{jj}=optimalcontroldynamicsleftside.term{jj};
                        end
                    end
                    varStruct.(varname).string{ii}{jj}=[varStruct.(varname).string{ii}{jj} ';'];
                else
                    varStruct.(varname).string{ii}{1}='out=[];';
                end
            end
            varStruct.(varname).multline=1;
            if arcnum.value>1
                varStruct.(varname).arcdependent=1;
            else
                varStruct.(varname).arcdependent=0;
                varStruct.(varname).string=varStruct.(varname).string{1};
            end
            varStruct.(varname).type='algterm';
        catch
            %lasterr
            varStruct.(varname).string=varname;
        end

    case 'OPTIMALCONTROLDYNAMICS'
        try
            arcfield=getarcclass(ocStruct,'canonicalsystemjacobian');
            arcnum=retrieveimpulsemodelinformation(ocStruct,'arcnum');
            varStruct.(varname).arcidentifier=field2arcidentifier(arcfield);
            for ii=1:numel(arcfield)
                optimalcontroldynamics=algebraicterm2string(ocStruct,'optimalcontroldynamics',1,arcfield{ii});
                if ~isempty(optimalcontroldynamics.term)
                    for jj=1:numel(optimalcontroldynamics.term)
                        if jj==1
                            varStruct.(varname).string{ii}{jj}=['dudt=' optimalcontroldynamics.term{jj}];
                        else
                            varStruct.(varname).string{ii}{jj}=optimalcontroldynamics.term{jj};
                        end
                    end
                    varStruct.(varname).string{ii}{jj}=[varStruct.(varname).string{ii}{jj} ';'];
                else
                    varStruct.(varname).string{ii}{1}='dudt=[];';
                end
            end
            varStruct.(varname).multline=1;
            varStruct.(varname).vectorize=1;
            if arcnum.value>1
                varStruct.(varname).arcdependent=3;
            else
                varStruct.(varname).arcdependent=0;
                varStruct.(varname).string=varStruct.(varname).string{1};
            end
            varStruct.(varname).type='algterm';
        catch
            %lasterr
            varStruct.(varname).string=varname;
        end

    case 'SYMBOLICOPTIMALCONTROLDYNAMICS'
        try
            arcfield=getarcclass(ocStruct,'canonicalsystemjacobian');
            arcnum=retrieveimpulsemodelinformation(ocStruct,'arcnum');
            varStruct.(varname).arcidentifier=field2arcidentifier(arcfield);
            for ii=1:numel(arcfield)
                optimalcontroldynamics=algebraicterm2string(ocStruct,'symbolicoptimalcontroldynamics',0,arcfield{ii});
                if ~isempty(optimalcontroldynamics.term)
                    for jj=1:numel(optimalcontroldynamics.term)
                        if jj==1
                            varStruct.(varname).string{ii}{jj}=['dudt=' optimalcontroldynamics.term{jj}];
                        else
                            varStruct.(varname).string{ii}{jj}=optimalcontroldynamics.term{jj};
                        end
                    end
                    varStruct.(varname).string{ii}{jj}=[varStruct.(varname).string{ii}{jj} ';'];
                else
                    varStruct.(varname).string{ii}{1}='dudt=[];';
                end
            end
            varStruct.(varname).multline=1;
            varStruct.(varname).vectorize=0;
            if arcnum.value>1
                varStruct.(varname).arcdependent=1;
            else
                varStruct.(varname).arcdependent=0;
                varStruct.(varname).string=varStruct.(varname).string{1};
            end
            varStruct.(varname).type='algterm';
        catch
            %lasterr
            varStruct.(varname).string=varname;
        end

    case 'TENSOROPTIMALCONTROLDYNAMICSLEFTSIDE'
        try
            arcfield=getarcclass(ocStruct,'canonicalsystemjacobian');
            arcnum=retrieveimpulsemodelinformation(ocStruct,'arcnum');
            varStruct.(varname).arcidentifier=field2arcidentifier(arcfield);
            for ii=1:numel(arcfield)
                optimalcontroldynamicsleftside=algebraicterm2string(ocStruct,'tensoroptimalcontroldynamicsleftside',0,arcfield{ii});
                if ~isempty(optimalcontroldynamicsleftside.term)
                    for jj=1:numel(optimalcontroldynamicsleftside.term)
                        if jj==1
                            varStruct.(varname).string{ii}{jj}=optimalcontroldynamicsleftside.term{jj};
                        else
                            varStruct.(varname).string{ii}{jj}=optimalcontroldynamicsleftside.term{jj};
                        end
                    end
                    %varStruct.(varname).string{ii}{jj}=[varStruct.(varname).string{ii}{jj} ';'];
                else
                    varStruct.(varname).string{ii}{1}='out=[];';
                end
            end
            varStruct.(varname).multline=1;
            if arcnum.value>1
                varStruct.(varname).arcdependent=1;
            else
                varStruct.(varname).arcdependent=0;
                varStruct.(varname).string=varStruct.(varname).string{1};
            end
            varStruct.(varname).type='algterm';
        catch
            %lasterr
            varStruct.(varname).string=varname;
        end

    case 'OPTIMALCONTROLDYNAMICSRIGHTSIDE'
        try
            pontryaginfunctionDuDX=algebraicterm2string(ocStruct,'pontryaginfunctionDuDX',0);
            for jj=1:numel(pontryaginfunctionDuDX.term)
                if jj==1
                    varStruct.(varname).string{jj}=pontryaginfunctionDuDX.term{jj};
                else
                    varStruct.(varname).string{jj}=pontryaginfunctionDuDX.term{jj};
                end
            end
            varStruct.(varname).string{jj}=[varStruct.(varname).string{jj} ';'];
            varStruct.(varname).multline=1;
            varStruct.(varname).arcdependent=0;
            varStruct.(varname).type='algterm';
        catch
            %lasterr
            varStruct.(varname).string=varname;
        end

    case 'CONSTRAINEDOPTIMALCONTROLDYNAMICSRIGHTSIDE'
        arcfield=getarcclass(ocStruct,'canonicalsystemjacobian');
        arcnum=retrieveimpulsemodelinformation(ocStruct,'arcnum');
        varStruct.(varname).arcidentifier=field2arcidentifier(arcfield);
        for ii=1:numel(arcfield)
            pontryaginfunctionDlmmcDX=algebraicterm2string(ocStruct,'pontryaginfunctionDlmmcDX',1,arcfield{ii});
            if ~isempty(pontryaginfunctionDlmmcDX.term)
                for jj=1:numel(pontryaginfunctionDlmmcDX.term)
                    if jj==1
                        varStruct.(varname).string{ii}{jj}=['tmp=' pontryaginfunctionDlmmcDX.term{jj}];
                    else
                        varStruct.(varname).string{ii}{jj}=pontryaginfunctionDlmmcDX.term{jj};
                    end
                end
                varStruct.(varname).string{ii}{jj}=[varStruct.(varname).string{ii}{jj} ';'];
            else
                varStruct.(varname).string{ii}{1}='tmp=[];';
            end
        end
        varStruct.(varname).multline=1;
        if arcnum.value>1
            varStruct.(varname).arcdependent=1;
        else
            varStruct.(varname).arcdependent=0;
            varStruct.(varname).string=varStruct.(varname).string{1};
        end
        varStruct.(varname).type='algterm';

    case 'JACOBIANOPTIMALCONTROLDYNAMICSRIGHTSIDE'
        try
            arcfield=getarcclass(ocStruct,'canonicalsystemjacobian');
            arcnum=retrieveimpulsemodelinformation(ocStruct,'arcnum');
            varStruct.(varname).arcidentifier=field2arcidentifier(arcfield);
            for ii=1:numel(arcfield)
                optimalcontroldynamicsleftside=algebraicterm2string(ocStruct,'jacobianoptimalcontroldynamicsrightside',0,arcfield{ii});
                if ~isempty(optimalcontroldynamicsleftside.term)
                    for jj=1:numel(optimalcontroldynamicsleftside.term)
                        if jj==1
                            varStruct.(varname).string{ii}{jj}=['out=' optimalcontroldynamicsleftside.term{jj}];
                        else
                            varStruct.(varname).string{ii}{jj}=optimalcontroldynamicsleftside.term{jj};
                        end
                    end
                    varStruct.(varname).string{ii}{jj}=[varStruct.(varname).string{ii}{jj} ';'];
                else
                    varStruct.(varname).string{ii}{1}='out=[];';
                end
            end
            varStruct.(varname).multline=1;
            if arcnum.value>1
                varStruct.(varname).arcdependent=1;
            else
                varStruct.(varname).arcdependent=0;
                varStruct.(varname).string=varStruct.(varname).string{1};
            end
            varStruct.(varname).type='algterm';
        catch
            %lasterr
            varStruct.(varname).string=varname;
        end

    case 'SYMBOLICCANONICALSYSTEMDYNAMICS'
        try
            canonicalsystem=algebraicterm2string(ocStruct,'symboliccanonicalsystem',0,'');
            for jj=1:numel(canonicalsystem.term)
                if jj==1
                    varStruct.(varname).string{jj}=['out=' canonicalsystem.term{jj}];
                else
                    varStruct.(varname).string{jj}=canonicalsystem.term{jj};
                end
            end
            varStruct.(varname).string{jj}=[varStruct.(varname).string{jj} ';'];
            varStruct.(varname).multline=1;
            varStruct.(varname).vectorize=0;
            varStruct.(varname).arcdependent=0;
            varStruct.(varname).type='algterm';
        catch
            %lasterr
            varStruct.(varname).string=varname;
        end

    case 'SYMBOLICEQUILIBRIUMEQUATION'
        try
            arcfield=getarcclass(ocStruct,'canonicalsystemjacobian');
            arcnum=retrieveimpulsemodelinformation(ocStruct,'arcnum');
            varStruct.(varname).arcidentifier=field2arcidentifier(arcfield);
            for ii=1:numel(arcfield)
                equilibriumequation=algebraicterm2string(ocStruct,'equilibriumequation',0,arcfield{ii});
                for jj=1:numel(equilibriumequation.term)
                    if jj==1
                        varStruct.(varname).string{ii}{jj}=['out=' equilibriumequation.term{jj}];
                    else
                        varStruct.(varname).string{ii}{jj}=equilibriumequation.term{jj};
                    end
                end
                varStruct.(varname).string{ii}{jj}=[varStruct.(varname).string{ii}{jj} ';'];
            end
            varStruct.(varname).multline=1;
            if arcnum.value>1
                varStruct.(varname).arcdependent=1;
            else
                varStruct.(varname).arcdependent=0;
                varStruct.(varname).string=varStruct.(varname).string{1};
            end
            varStruct.(varname).type='algterm';
        catch
            %lasterr
            varStruct.(varname).string=varname;
        end

    case 'SYMBOLICCANONICALSYSTEMALGEBRAICEQUATION'
        try
            arcfield=getarcclass(ocStruct,'algebraicequation');
            varStruct.(varname).arcidentifier=field2arcidentifier(arcfield);
            arcnum=retrieveimpulsemodelinformation(ocStruct,'arcnum');
            for ii=1:numel(arcfield)
                algebraicequation=algebraicterm2string(ocStruct,'symbolicalgebraicequation',0,arcfield{ii});
                if ~isempty(algebraicequation.term)
                    for jj=1:numel(algebraicequation.term)
                        if jj==1
                            varStruct.(varname).string{ii}{jj}=['outae=' algebraicequation.term{jj}];
                        else
                            varStruct.(varname).string{ii}{jj}=algebraicequation.term{jj};
                        end
                    end
                    varStruct.(varname).string{ii}{jj}=[varStruct.(varname).string{ii}{jj} ';'];
                else
                    varStruct.(varname).string{ii}{1}='outae=sym([]);';
                end
            end
            varStruct.(varname).vectorize=0;
            varStruct.(varname).multline=1;
            if arcnum.value>1
                varStruct.(varname).arcdependent=1;
            else
                varStruct.(varname).arcdependent=0;
                varStruct.(varname).string=varStruct.(varname).string{1};
            end
            varStruct.(varname).type='algterm';
        catch
            %lasterr
            varStruct.(varname).string=varname;
        end

    case 'CANONICALSYSTEMJACOBIAN'
        try
            arcfield=getarcclass(ocStruct,'canonicalsystemjacobian');
            arcnum=retrieveimpulsemodelinformation(ocStruct,'arcnum');
            varStruct.(varname).arcidentifier=field2arcidentifier(arcfield);
            for ii=1:numel(arcfield)
                canonicalsystemjacobian=algebraicterm2string(ocStruct,'canonicalsystemjacobian',0,arcfield{ii});
                for jj=1:numel(canonicalsystemjacobian.term)
                    if jj==1
                        varStruct.(varname).string{ii}{jj}=['out=' canonicalsystemjacobian.term{jj}];
                    else
                        varStruct.(varname).string{ii}{jj}=canonicalsystemjacobian.term{jj};
                    end
                end
                varStruct.(varname).string{ii}{jj}=[varStruct.(varname).string{ii}{jj} ';'];
            end
            varStruct.(varname).multline=1;
            if arcnum.value>1
                varStruct.(varname).arcdependent=1;
            else
                varStruct.(varname).arcdependent=0;
                varStruct.(varname).string=varStruct.(varname).string{1};
            end
            varStruct.(varname).type='algterm';
        catch
            %lasterr
            varStruct.(varname).string=varname;
        end

    case 'STATECONTROLCANONICALSYSTEMJACOBIAN'
        try
            arcfield=getarcclass(ocStruct,'canonicalsystemjacobian');
            arcnum=retrieveimpulsemodelinformation(ocStruct,'arcnum');
            varStruct.(varname).arcidentifier=field2arcidentifier(arcfield);
            for ii=1:numel(arcfield)
                canonicalsystemjacobian=algebraicterm2string(ocStruct,'canonicalsystemjacobianstatecontrol',0,arcfield{ii});
                for jj=1:numel(canonicalsystemjacobian.term)
                    if jj==1
                        varStruct.(varname).string{ii}{jj}=['out=' canonicalsystemjacobian.term{jj}];
                    else
                        varStruct.(varname).string{ii}{jj}=canonicalsystemjacobian.term{jj};
                    end
                end
                varStruct.(varname).string{ii}{jj}=[varStruct.(varname).string{ii}{jj} ';'];
            end
            varStruct.(varname).multline=1;
            if arcnum.value>1
                varStruct.(varname).arcdependent=1;
            else
                varStruct.(varname).arcdependent=0;
                varStruct.(varname).string=varStruct.(varname).string{1};
            end
            varStruct.(varname).type='algterm';
        catch
            %lasterr
            varStruct.(varname).string=varname;
        end

    case 'STATECOSTATEJACOBIAN'
        try
            arcfield=getarcclass(ocStruct,'canonicalsystemjacobian');
            arcnum=retrieveimpulsemodelinformation(ocStruct,'arcnum');
            varStruct.(varname).arcidentifier=field2arcidentifier(arcfield);
            for ii=1:numel(arcfield)
                statecostatejacobian=algebraicterm2string(ocStruct,'statecostatejacobian',0,arcfield{ii});
                for jj=1:numel(statecostatejacobian.term)
                    if jj==1
                        varStruct.(varname).string{ii}{jj}=['out=' statecostatejacobian.term{jj}];
                    else
                        varStruct.(varname).string{ii}{jj}=statecostatejacobian.term{jj};
                    end
                end
                varStruct.(varname).string{ii}{jj}=[varStruct.(varname).string{ii}{jj} ';'];
            end
            varStruct.(varname).multline=1;
            if arcnum.value>1
                varStruct.(varname).arcdependent=1;
            else
                varStruct.(varname).arcdependent=0;
                varStruct.(varname).string=varStruct.(varname).string{1};
            end
            varStruct.(varname).type='algterm';
        catch
            %lasterr
            varStruct.(varname).string=varname;
        end

    case 'SPECIFICCANONICALSYSTEMDYNAMICS'
        try
            arcfield=getarcclass(ocStruct,'canonicalsystemjacobian');
            arcnum=retrieveimpulsemodelinformation(ocStruct,'arcnum');
            varStruct.(varname).arcidentifier=field2arcidentifier(arcfield);
            for ii=1:numel(arcfield)
                statecostatejacobian=algebraicterm2string(ocStruct,'specificcanonicalsystem',1,arcfield{ii});
                for jj=1:numel(statecostatejacobian.term)
                    if jj==1
                        varStruct.(varname).string{ii}{jj}=['out=' statecostatejacobian.term{jj}];
                    else
                        varStruct.(varname).string{ii}{jj}=statecostatejacobian.term{jj};
                    end
                end
                varStruct.(varname).string{ii}{jj}=[varStruct.(varname).string{ii}{jj} ';'];
            end
            varStruct.(varname).multline=1;
            if arcnum.value>1
                varStruct.(varname).arcdependent=1;
            else
                varStruct.(varname).arcdependent=0;
                varStruct.(varname).string=varStruct.(varname).string{1};
            end
            varStruct.(varname).type='algterm';
        catch
            %lasterr
            varStruct.(varname).string=varname;
        end

    case 'SPECIFICPONTRYAGINFUNCTION'
        try
            arcfield=getarcclass(ocStruct,'canonicalsystemjacobian');
            arcnum=retrieveimpulsemodelinformation(ocStruct,'arcnum');
            varStruct.(varname).arcidentifier=field2arcidentifier(arcfield);
            for ii=1:numel(arcfield)
                specificpontryaginfunction=algebraicterm2string(ocStruct,'specificpontryaginfunction',1,arcfield{ii});
                varStruct.(varname).string{ii}=['out=' specificpontryaginfunction.term ';'];
            end
            varStruct.(varname).multline=0;
            if 1%arcnum.value>1
                varStruct.(varname).arcdependent=1;
            else
                varStruct.(varname).arcdependent=0;
                varStruct.(varname).string=varStruct.(varname).string{1};
            end
            varStruct.(varname).type='algterm';
        catch
            %lasterr
            varStruct.(varname).string=varname;
        end

    case 'SYMBOLICCANONICALSYSTEMJACOBIAN'
        try
            arcfield=getarcclass(ocStruct,'canonicalsystemjacobian');
            arcnum=retrieveimpulsemodelinformation(ocStruct,'arcnum');
            varStruct.(varname).arcidentifier=field2arcidentifier(arcfield);
            for ii=1:numel(arcfield)
                canonicalsystemjacobian=algebraicterm2string(ocStruct,'symboliccanonicalsystemjacobian',0,arcfield{ii});
                for jj=1:numel(canonicalsystemjacobian.term)
                    if jj==1
                        varStruct.(varname).string{ii}{jj}=['out=' canonicalsystemjacobian.term{jj}];
                    else
                        varStruct.(varname).string{ii}{jj}=canonicalsystemjacobian.term{jj};
                    end
                end
                varStruct.(varname).string{ii}{jj}=[varStruct.(varname).string{ii}{jj} ';'];
            end
            varStruct.(varname).multline=1;
            if arcnum.value>1
                varStruct.(varname).arcdependent=1;
            else
                varStruct.(varname).arcdependent=0;
                varStruct.(varname).string=varStruct.(varname).string{1};
            end
            varStruct.(varname).type='algterm';
        catch
            %lasterr
            varStruct.(varname).string=varname;
        end

    case 'SYMBOLICSTATECONTROLCANONICALSYSTEMJACOBIAN'
        try
            arcfield=getarcclass(ocStruct,'canonicalsystemjacobian');
            arcnum=retrieveimpulsemodelinformation(ocStruct,'arcnum');
            varStruct.(varname).arcidentifier=field2arcidentifier(arcfield);
            for ii=1:numel(arcfield)
                canonicalsystemjacobian=algebraicterm2string(ocStruct,'symbolicstatecontrolcanonicalsystemjacobian',0,arcfield{ii});
                for jj=1:numel(canonicalsystemjacobian.term)
                    if jj==1
                        varStruct.(varname).string{ii}{jj}=['out=' canonicalsystemjacobian.term{jj}];
                    else
                        varStruct.(varname).string{ii}{jj}=canonicalsystemjacobian.term{jj};
                    end
                end
                varStruct.(varname).string{ii}{jj}=[varStruct.(varname).string{ii}{jj} ';'];
            end
            varStruct.(varname).multline=1;
            if arcnum.value>1
                varStruct.(varname).arcdependent=1;
            else
                varStruct.(varname).arcdependent=0;
                varStruct.(varname).string=varStruct.(varname).string{1};
            end
            varStruct.(varname).type='algterm';
        catch
            %lasterr
            varStruct.(varname).string=varname;
        end


    case 'SYMBOLICSTATECONTROLCANONICALSYSTEMDYNAMICS'
         try
            arcfield=getarcclass(ocStruct,'canonicalsystemjacobian');
            arcnum=retrieveimpulsemodelinformation(ocStruct,'arcnum');
            varStruct.(varname).arcidentifier=field2arcidentifier(arcfield);
            for ii=1:numel(arcfield)
                canonicalsystemjacobian=algebraicterm2string(ocStruct,'symbolicstatecontrolcanonicalsystem',0,arcfield{ii});
                for jj=1:numel(canonicalsystemjacobian.term)
                    if jj==1
                        varStruct.(varname).string{ii}{jj}=['out=' canonicalsystemjacobian.term{jj}];
                    else
                        varStruct.(varname).string{ii}{jj}=canonicalsystemjacobian.term{jj};
                    end
                end
                varStruct.(varname).string{ii}{jj}=[varStruct.(varname).string{ii}{jj} ';'];
            end
            varStruct.(varname).multline=1;
            if arcnum.value>1
                varStruct.(varname).arcdependent=1;
            else
                varStruct.(varname).arcdependent=0;
                varStruct.(varname).string=varStruct.(varname).string{1};
            end
            varStruct.(varname).type='algterm';
        catch
            %lasterr
            varStruct.(varname).string=varname;
        end

    case 'CANONICALSYSTEMPARAMETERJACOBIAN'
        try
            arcfield=getarcclass(ocStruct,'canonicalsystemjacobian');
            varStruct.(varname).arcidentifier=field2arcidentifier(arcfield);
            for ii=1:numel(arcfield)
                canonicalsystemparameterjacobian=algebraicterm2string(ocStruct,'canonicalsystemparameterjacobian',0,arcfield{ii});
                for jj=1:numel(canonicalsystemparameterjacobian.term)
                    if jj==1
                        varStruct.(varname).string{ii}{jj}=['out=' canonicalsystemparameterjacobian.term{jj}];
                    else
                        varStruct.(varname).string{ii}{jj}=canonicalsystemparameterjacobian.term{jj};
                    end
                end
                varStruct.(varname).string{ii}{jj}=[varStruct.(varname).string{ii}{jj} ';'];
            end
            varStruct.(varname).multline=1;
            arcnum=retrieveimpulsemodelinformation(ocStruct,'arcnum');
            if arcnum.value>1
                varStruct.(varname).arcdependent=1;
            else
                varStruct.(varname).arcdependent=0;
                varStruct.(varname).string=varStruct.(varname).string{1};
            end
            varStruct.(varname).type='algterm';
        catch
            %lasterr
            varStruct.(varname).string=varname;
        end

    case 'STATECONTROLCANONICALSYSTEMPARAMETERJACOBIAN'
        try
            arcfield=getarcclass(ocStruct,'canonicalsystemjacobian');
            varStruct.(varname).arcidentifier=field2arcidentifier(arcfield);
            for ii=1:numel(arcfield)
                canonicalsystemparameterjacobian=algebraicterm2string(ocStruct,'statecontrolcanonicalsystemparameterjacobian',0,arcfield{ii});
                for jj=1:numel(canonicalsystemparameterjacobian.term)
                    if jj==1
                        varStruct.(varname).string{ii}{jj}=['out=' canonicalsystemparameterjacobian.term{jj}];
                    else
                        varStruct.(varname).string{ii}{jj}=canonicalsystemparameterjacobian.term{jj};
                    end
                end
                varStruct.(varname).string{ii}{jj}=[varStruct.(varname).string{ii}{jj} ';'];
            end
            varStruct.(varname).multline=1;
            arcnum=retrieveimpulsemodelinformation(ocStruct,'arcnum');
            if arcnum.value>1
                varStruct.(varname).arcdependent=1;
            else
                varStruct.(varname).arcdependent=0;
                varStruct.(varname).string=varStruct.(varname).string{1};
            end
            varStruct.(varname).type='algterm';
        catch
            %lasterr
            varStruct.(varname).string=varname;
        end

    case 'CANONICALSYSTEMDERIVATIVETIME'
        try
            arcfield=getarcclass(ocStruct,'canonicalsystemjacobian');
            varStruct.(varname).arcidentifier=field2arcidentifier(arcfield);
            for ii=1:numel(arcfield)
                canonicalsystemderivativetime=algebraicterm2string(ocStruct,'canonicalsystemderivativetime',0,arcfield{ii});
                for jj=1:numel(canonicalsystemderivativetime.term)
                    if jj==1
                        varStruct.(varname).string{ii}{jj}=['out=' canonicalsystemderivativetime.term{jj}];
                    else
                        varStruct.(varname).string{ii}{jj}=canonicalsystemderivativetime.term{jj};
                    end
                end
                varStruct.(varname).string{ii}{jj}=[varStruct.(varname).string{ii}{jj} ';'];
            end
            varStruct.(varname).multline=1;
            arcnum=retrieveimpulsemodelinformation(ocStruct,'arcnum');
            if arcnum.value>1
                varStruct.(varname).arcdependent=1;
            else
                varStruct.(varname).arcdependent=0;
                varStruct.(varname).string=varStruct.(varname).string{1};
            end
            varStruct.(varname).type='algterm';
        catch
            %lasterr
            varStruct.(varname).string=varname;
        end

    case 'ARCIDENTIFIERALGEBRAICEQUATION'
        varStruct.(varname).arcdependent=1;
        try
            arcfield=getarcclass(ocStruct,'algebraicequation');
            for ii=1:numel(arcfield)
                arcidentifier=field2arcidentifier(arcfield{ii});
                if numel(arcidentifier)>1
                    varStruct.(varname).string{ii}=['{' arcidentifier{1}];
                    for jj=2:numel(arcidentifier)
                        varStruct.(varname).string{ii}=[varStruct.(varname).string{ii} ',' arcidentifier{jj}];
                    end
                    varStruct.(varname).string{ii}=[varStruct.(varname).string{ii} '}'];
                else
                    varStruct.(varname).string{ii}=arcidentifier{1};
                end
            end
        catch
            %lasterr
            varStruct.(varname).string={'0'};
        end

    case 'ARCIDENTIFIERCANONICALSYSTEMJACOBIAN'
        varStruct.(varname).arcdependent=1;
        try
            arcfield=getarcclass(ocStruct,'canonicalsystemjacobian');
            for ii=1:numel(arcfield)
                arcidentifier=field2arcidentifier(arcfield{ii});
                if numel(arcidentifier)>1
                    varStruct.(varname).string{ii}=['{' arcidentifier{1}];
                    for jj=2:numel(arcidentifier)
                        varStruct.(varname).string{ii}=[varStruct.(varname).string{ii} ',' arcidentifier{jj}];
                    end
                    varStruct.(varname).string{ii}=[varStruct.(varname).string{ii} '}'];
                else
                    varStruct.(varname).string{ii}=arcidentifier{1};
                end
            end
        catch
            %lasterr
            varStruct.(varname).string={'0'};
        end

    case 'CANONICALSYSTEMHESSIAN'
        try
            arcfield=getarcclass(ocStruct,'canonicalsystemjacobian');
            varStruct.(varname).arcidentifier=field2arcidentifier(arcfield);
            for ii=1:numel(arcfield)
                canonicalsystemhessian=algebraicterm2string(ocStruct,'canonicalsystemhessian',0,arcfield{ii});
                numlines=canonicalsystemhessian.info;
                counter=0;
                for jj=1:numel(canonicalsystemhessian.term)
                    if rem(jj,numlines)==1
                        counter=counter+1;
                        varStruct.(varname).string{ii}{jj}=['out(:,:,' num2str(counter) ')=' canonicalsystemhessian.term{jj}];
                    elseif rem(jj,numlines)==0
                        varStruct.(varname).string{ii}{jj}=[canonicalsystemhessian.term{jj} ';'];
                    else
                        varStruct.(varname).string{ii}{jj}=canonicalsystemhessian.term{jj};
                    end
                end
            end
            varStruct.(varname).multline=1;
            arcnum=retrieveimpulsemodelinformation(ocStruct,'arcnum');
            if arcnum.value>1
                varStruct.(varname).arcdependent=1;
            else
                varStruct.(varname).arcdependent=0;
                varStruct.(varname).string=varStruct.(varname).string{1};
            end
            varStruct.(varname).type='algterm';
        catch
            %lasterr
            varStruct.(varname).string=varname;
        end

    case 'CANONICALSYSTEMTOTALHESSIAN'
        try
            arcfield=getarcclass(ocStruct,'canonicalsystemjacobian');
            varStruct.(varname).arcidentifier=field2arcidentifier(arcfield);
            for ii=1:numel(arcfield)
                canonicalsystemparameterhessian=algebraicterm2string(ocStruct,'canonicalsystemtotalhessian',0,arcfield{ii});
                numlines=canonicalsystemparameterhessian.info;
                counter=0;
                for jj=1:numel(canonicalsystemparameterhessian.term)
                    if rem(jj,numlines)==1
                        counter=counter+1;
                        varStruct.(varname).string{ii}{jj}=['out(:,:,' num2str(counter) ')=' canonicalsystemparameterhessian.term{jj}];
                    elseif rem(jj,numlines)==0
                        varStruct.(varname).string{ii}{jj}=[canonicalsystemparameterhessian.term{jj} ';'];
                    else
                        varStruct.(varname).string{ii}{jj}=canonicalsystemparameterhessian.term{jj};
                    end
                end
            end
            varStruct.(varname).multline=1;
            arcnum=retrieveimpulsemodelinformation(ocStruct,'arcnum');
            if arcnum.value>1
                varStruct.(varname).arcdependent=1;
            else
                varStruct.(varname).arcdependent=0;
                varStruct.(varname).string=varStruct.(varname).string{1};
            end
            varStruct.(varname).type='algterm';
        catch
            %lasterr
            varStruct.(varname).string=varname;
        end

    case 'D2HAMILTONIANDU2'
        try
            arcfield=getarcclass(ocStruct,'canonicalsystemjacobian');
            varStruct.(varname).arcidentifier=field2arcidentifier(arcfield);
            for ii=1:numel(arcfield)
                pontryaginfunctionDu2=algebraicterm2string(ocStruct,'pontryaginfunctionDu2',0,arcfield{ii});
                counter=0;
                for jj=1:numel(pontryaginfunctionDu2.term)
                    if jj==1
                        counter=counter+1;
                        varStruct.(varname).string{ii}{jj}=['out=' pontryaginfunctionDu2.term{jj}];
                    else
                        varStruct.(varname).string{ii}{jj}=pontryaginfunctionDu2.term{jj};
                    end
                end
                varStruct.(varname).string{ii}{jj}=[varStruct.(varname).string{ii}{jj} ';'];
            end
            varStruct.(varname).multline=1;
            arcnum=retrieveimpulsemodelinformation(ocStruct,'arcnum');
            if arcnum.value>1
                varStruct.(varname).arcdependent=1;
            else
                varStruct.(varname).arcdependent=0;
                varStruct.(varname).string=varStruct.(varname).string{1};
            end
            varStruct.(varname).type='algterm';
        catch
            %lasterr
            varStruct.(varname).string=varname;
        end

    case 'PARAMETERHESSIANCOORD'
        arcfield=getarcclass(ocStruct,'canonicalsystemjacobian');
        varStruct.(varname).arcidentifier=field2arcidentifier(arcfield);
        parameternum=retrieveimpulsemodelinformation(ocStruct,'parameternum');
        for ii=1:numel(arcfield)
            equationnum=retrieveimpulsemodelinformation(ocStruct,'canonicalsystemequationnum',field2arcidentifier(arcfield{ii}));
            varStruct.(varname).string{ii}=['1:' int2str(equationnum.value) ',1:'  int2str(equationnum.value) ',' int2str(equationnum.value+1) ':' int2str(parameternum.value+equationnum.value)];
            arcnum=retrieveimpulsemodelinformation(ocStruct,'arcnum');
            if arcnum.value>1
                varStruct.(varname).arcdependent=1;
            else
                varStruct.(varname).arcdependent=0;
                varStruct.(varname).string=varStruct.(varname).string{1};
            end
        end

    case 'ARCIDENTIFIERLAGRANGEMULTIPLIER'
        varStruct.(varname).arcdependent=1;
        try
            arcfield=getarcclass(ocStruct,'lagrangemultiplier');
            for ii=1:numel(arcfield)
                arcidentifier=field2arcidentifier(arcfield{ii});
                if numel(arcidentifier)>1
                    varStruct.(varname).string{ii}=['{' arcidentifier{1}];
                    for jj=2:numel(arcidentifier)
                        varStruct.(varname).string{ii}=[varStruct.(varname).string{ii} ',' arcidentifier{jj}];
                    end
                    varStruct.(varname).string{ii}=[varStruct.(varname).string{ii} '}'];
                else
                    varStruct.(varname).string{ii}=arcidentifier{1};
                end
            end
        catch
            %lasterr
            varStruct.(varname).string={'0'};
        end

    case 'OPTIMALCONTROLVALUE'
        try
            arcfield=getarcclass(ocStruct,'controlvalue');
            varStruct.(varname).arcidentifier=field2arcidentifier(arcfield);
            for ii=1:numel(arcfield)
                controlvalue=algebraicterm2string(ocStruct,'controlvalue',1,arcfield{ii});
                if ~isempty(controlvalue.term)
                    for jj=1:numel(controlvalue.term)
                        if jj==1
                            varStruct.(varname).string{ii}{jj}=['out=' controlvalue.term{jj}];
                        else
                            varStruct.(varname).string{ii}{jj}=controlvalue.term{jj};
                        end
                    end
                    varStruct.(varname).string{ii}{jj}=[varStruct.(varname).string{ii}{jj} ';'];
                else
                    varStruct.(varname).string{ii}{1}='out=[];';
                end
            end
            varStruct.(varname).multline=1;
            varStruct.(varname).arcdependent=1;
            varStruct.(varname).vectorize=1;
            varStruct.(varname).type='algterm';
        catch
            %lasterr
            varStruct.(varname).string=varname;
        end

    case 'OPTIMALIMPULSECONTROLVALUE'
        try
            arcfield=getarcclass(ocStruct,'impulsecontrolvalue');
            arcfield=arcfield(1);
            varStruct.(varname).arcidentifier=field2arcidentifier(arcfield);
            for ii=1:numel(arcfield)
                controlvalue=algebraicterm2string(ocStruct,'impulsecontrolvalue',1,arcfield{ii});
                for jj=1:numel(controlvalue.term)
                    if jj==1
                        varStruct.(varname).string{ii}{jj}=['out=' controlvalue.term{jj}];
                    else
                        varStruct.(varname).string{ii}{jj}=controlvalue.term{jj};
                    end
                end
                varStruct.(varname).string{ii}{jj}=[varStruct.(varname).string{ii}{jj} ';'];
            end
            varStruct.(varname).multline=1;
            varStruct.(varname).vectorize=0;
            varStruct.(varname).type='algterm';
            if numel(arcfield)>1
                varStruct.(varname).arcdependent=1;
            else
                varStruct.(varname).arcdependent=0;
                varStruct.(varname).string=varStruct.(varname).string{1};
            end
        catch
            %lasterr
            varStruct.(varname).string=varname;
        end

    case 'SYMBOLICOPTIMALIMPULSECONTROLVALUE'
        try
            arcfield=getarcclass(ocStruct,'impulsecontrolvalue');
            arcfield=arcfield(1);
            varStruct.(varname).arcidentifier=field2arcidentifier(arcfield);
            for ii=1:numel(arcfield)
                symbolicimpulsecontrolvalue=algebraicterm2string(ocStruct,'symbolicimpulsecontrolvalue',0,arcfield{ii});
                for jj=1:numel(symbolicimpulsecontrolvalue.term)
                    if jj==1
                        varStruct.(varname).string{ii}{jj}=['out=' symbolicimpulsecontrolvalue.term{jj}];
                    else
                        varStruct.(varname).string{ii}{jj}=symbolicimpulsecontrolvalue.term{jj};
                    end
                end
                varStruct.(varname).string{ii}{jj}=[varStruct.(varname).string{ii}{jj} ';'];
            end
            varStruct.(varname).multline=1;
            varStruct.(varname).vectorize=0;
            varStruct.(varname).type='algterm';
            if numel(arcfield)>1
                varStruct.(varname).arcdependent=1;
            else
                varStruct.(varname).arcdependent=0;
                varStruct.(varname).string=varStruct.(varname).string{1};
            end
        catch
            %lasterr
            varStruct.(varname).string=varname;
        end
    case 'OPTIMALCONTROLVALUEIMPLICIT'
        try
            arcfield=getarcclass(ocStruct,'controlvalue');
            varStruct.(varname).arcidentifier=field2arcidentifier(arcfield);
            for ii=1:numel(arcfield)
                controlvalue=algebraicterm2string(ocStruct,'controlvalue',1,arcfield{ii});
                for jj=1:numel(controlvalue.term)
                    if jj==1
                        varStruct.(varname).string{ii}{jj}=['out=' controlvalue.term{jj}];
                    else
                        varStruct.(varname).string{ii}{jj}=controlvalue.term{jj};
                    end
                end
                varStruct.(varname).string{ii}{jj}=[varStruct.(varname).string{ii}{jj} ';'];
            end
            varStruct.(varname).multline=1;
            varStruct.(varname).arcdependent=2;
            varStruct.(varname).vectorize=2;
            varStruct.(varname).type='algterm';
        catch
            %lasterr
            varStruct.(varname).string=varname;
        end

    case 'OPTIMALCONTROLVALUE4STATECONTROL'
        try
            arcfield=getarcclass(ocStruct,'controlvalue');
            varStruct.(varname).arcidentifier=field2arcidentifier(arcfield);
            for ii=1:numel(arcfield)
                controlvalue=algebraicterm2string(ocStruct,'optimalvalue4statecontrol',1,arcfield{ii});
                for jj=1:numel(controlvalue.term)
                    if jj==1
                        varStruct.(varname).string{ii}{jj}=['out=' controlvalue.term{jj}];
                    else
                        varStruct.(varname).string{ii}{jj}=controlvalue.term{jj};
                    end
                end
                varStruct.(varname).string{ii}{jj}=[varStruct.(varname).string{ii}{jj} ';'];
            end
            varStruct.(varname).multline=1;
            varStruct.(varname).arcdependent=1;
            varStruct.(varname).vectorize=1;
            varStruct.(varname).type='algterm';
        catch
            %lasterr
            varStruct.(varname).string=varname;
        end

    case 'OPTIMALCOSTATEVALUE'
        try
            arcfield=getarcclass(ocStruct,'controlvalue');
            arcnum=retrieveimpulsemodelinformation(ocStruct,'arcnum');
            varStruct.(varname).arcidentifier=field2arcidentifier(arcfield);
            for ii=1:numel(arcfield)
                optimalcostatevalue=algebraicterm2string(ocStruct,'optimalcostatevalue',1,arcfield{ii});
                for jj=1:numel(optimalcostatevalue.term)
                    if jj==1
                        varStruct.(varname).string{ii}{jj}=['out=' optimalcostatevalue.term{jj}];
                    else
                        varStruct.(varname).string{ii}{jj}=optimalcostatevalue.term{jj};
                    end
                end
                varStruct.(varname).string{ii}{jj}=[varStruct.(varname).string{ii}{jj} ';'];
            end
            varStruct.(varname).multline=1;
            if arcnum.value>1
                varStruct.(varname).arcdependent=1;
            else
                varStruct.(varname).arcdependent=0;
                varStruct.(varname).string=varStruct.(varname).string{1};
            end
            varStruct.(varname).vectorize=1;
            varStruct.(varname).type='algterm';
        catch
            %lasterr
            varStruct.(varname).string=varname;
        end

    case 'EXPLICITSTATEVALUE'
        try
            arcfield=getarcclass(ocStruct,'controlvalue');
            arcnum=retrieveimpulsemodelinformation(ocStruct,'arcnum');
            varStruct.(varname).arcidentifier=field2arcidentifier(arcfield);
            for ii=1:numel(arcfield)
                explicitstatevalue=algebraicterm2string(ocStruct,'explicitstatevalue',1,arcfield{ii});
                for jj=1:numel(explicitstatevalue.term)
                    if jj==1
                        varStruct.(varname).string{ii}{jj}=['out=' explicitstatevalue.term{jj}];
                    else
                        varStruct.(varname).string{ii}{jj}=explicitstatevalue.term{jj};
                    end
                end
                varStruct.(varname).string{ii}{jj}=[varStruct.(varname).string{ii}{jj} ';'];
            end
            varStruct.(varname).multline=1;
            if arcnum.value>1
                varStruct.(varname).arcdependent=1;
            else
                varStruct.(varname).arcdependent=0;
                varStruct.(varname).string=varStruct.(varname).string{1};
            end
            varStruct.(varname).vectorize=1;
            varStruct.(varname).type='algterm';
        catch
            %lasterr
            varStruct.(varname).string=varname;
        end

    case 'SYMBOLICSTATEVALUE'
        try
            arcfield=getarcclass(ocStruct,'controlvalue');
            varStruct.(varname).arcidentifier=field2arcidentifier(arcfield);
            for ii=1:numel(arcfield)
                symbolicstatevalue=algebraicterm2string(ocStruct,'symbolicstatevalue',0,arcfield{ii});
                for jj=1:numel(symbolicstatevalue.term)
                    if jj==1
                        varStruct.(varname).string{ii}{jj}=['out=' symbolicstatevalue.term{jj}];
                    else
                        varStruct.(varname).string{ii}{jj}=symbolicstatevalue.term{jj};
                    end
                end
                varStruct.(varname).string{ii}{jj}=[varStruct.(varname).string{ii}{jj} ';'];
            end
            varStruct.(varname).multline=1;
            varStruct.(varname).arcdependent=1;
            varStruct.(varname).vectorize=0;
            varStruct.(varname).type='algterm';
        catch
            %lasterr
            varStruct.(varname).string=varname;
        end

    case 'SYMBOLICOPTIMALCONTROLVALUE'
        try
            arcfield=getarcclass(ocStruct,'controlvalue');
            varStruct.(varname).arcidentifier=field2arcidentifier(arcfield);
            for ii=1:numel(arcfield)
                controlvalue=algebraicterm2string(ocStruct,'symboliccontrolvalue',0,arcfield{ii});
                if ~isempty(controlvalue.term)
                    for jj=1:numel(controlvalue.term)
                        if jj==1
                            varStruct.(varname).string{ii}{jj}=['out=' controlvalue.term{jj}];
                        else
                            varStruct.(varname).string{ii}{jj}=controlvalue.term{jj};
                        end
                    end
                    varStruct.(varname).string{ii}{jj}=[varStruct.(varname).string{ii}{jj} ';'];
                else
                    varStruct.(varname).string{ii}='out=sym([]);';
                end
            end
            if ~isempty(controlvalue.term)
                varStruct.(varname).multline=1;
            else
                varStruct.(varname).multline=0;
            end
            varStruct.(varname).arcdependent=1;
            varStruct.(varname).vectorize=0;
            varStruct.(varname).type='algterm';
        catch
            %lasterr
            varStruct.(varname).string=varname;
        end

    case 'SYMBOLICOPTIMALCONTROLVALUE4STATECONTROL'
        try
            arcfield=getarcclass(ocStruct,'controlvalue');
            varStruct.(varname).arcidentifier=field2arcidentifier(arcfield);
            for ii=1:numel(arcfield)
                controlvalue=algebraicterm2string(ocStruct,'symboliccontrolvalue4statecontrol',0,arcfield{ii});
                for jj=1:numel(controlvalue.term)
                    if jj==1
                        varStruct.(varname).string{ii}{jj}=['out=' controlvalue.term{jj}];
                    else
                        varStruct.(varname).string{ii}{jj}=controlvalue.term{jj};
                    end
                end
                varStruct.(varname).string{ii}{jj}=[varStruct.(varname).string{ii}{jj} ';'];
            end
            varStruct.(varname).multline=1;
            varStruct.(varname).arcdependent=1;
            varStruct.(varname).vectorize=0;
            varStruct.(varname).type='algterm';
        catch
            %lasterr
            varStruct.(varname).string=varname;
        end

    case 'SYMBOLICOPTIMALCOSTATEVALUE'
        try
            arcfield=getarcclass(ocStruct,'controlvalue');
            varStruct.(varname).arcidentifier=field2arcidentifier(arcfield);
            for ii=1:numel(arcfield)
                controlvalue=algebraicterm2string(ocStruct,'symboliccostatevalue',0,arcfield{ii});
                for jj=1:numel(controlvalue.term)
                    if jj==1
                        varStruct.(varname).string{ii}{jj}=['out=' controlvalue.term{jj}];
                    else
                        varStruct.(varname).string{ii}{jj}=controlvalue.term{jj};
                    end
                end
                varStruct.(varname).string{ii}{jj}=[varStruct.(varname).string{ii}{jj} ';'];
            end
            varStruct.(varname).multline=1;
            varStruct.(varname).arcdependent=1;
            varStruct.(varname).vectorize=0;
            varStruct.(varname).type='algterm';
        catch
            %lasterr
            varStruct.(varname).string=varname;
        end

    case 'CONTROLLAGRANGEMULTIPLIER'
        try
            arcid=getarcclass(ocStruct,'lagrangemultcc');
            varStruct.(varname).arcidentifier=field2arcidentifier(arcid);
            for ii=1:numel(arcid)
                lagrangemultcc=algebraicterm2string(ocStruct,'lagrangemultcc',1,arcid{ii});
                for jj=1:numel(lagrangemultcc.term)
                    if jj==1
                        varStruct.(varname).string{ii}{jj}=['out=' lagrangemultcc.term{jj}];
                    else
                        varStruct.(varname).string{ii}{jj}=lagrangemultcc.term{jj};
                    end
                end
                varStruct.(varname).string{ii}{jj}=[varStruct.(varname).string{ii}{jj} ';'];
            end
            varStruct.(varname).arcdependent=1;
            varStruct.(varname).multline=1;
            varStruct.(varname).vectorize=1;
            varStruct.(varname).type='algterm';
        catch
            %lasterr
            varStruct.(varname).string=varname;
        end

    case 'STATELAGRANGEMULTIPLIER'
        try
            arcid=getarcclass(ocStruct,'lagrangemultsc');
            varStruct.(varname).arcidentifier=field2arcidentifier(arcid);
            if isempty(arcid)
                varStruct.(varname).string='';
            else
                for ii=1:numel(arcid)
                    lagrangemultsc=algebraicterm2string(ocStruct,'lagrangemultsc',1,arcid{ii});
                    for jj=1:numel(lagrangemultsc.term)
                        if jj==1
                            varStruct.(varname).string{ii}{jj}=['out=' lagrangemultsc.term{jj}];
                        else
                            varStruct.(varname).string{ii}{jj}=lagrangemultsc.term{jj};
                        end
                    end
                    varStruct.(varname).string{ii}{jj}=[varStruct.(varname).string{ii}{jj} ';'];
                end
            end
            varStruct.(varname).arcdependent=1;
            varStruct.(varname).multline=1;
            varStruct.(varname).vectorize=1;
            varStruct.(varname).type='algterm';
        catch
            %%lasterr
            varStruct.(varname).string=varname;
        end

    case 'SYMBOLICCONTROLLAGRANGEMULTIPLIER'
        try
            arcid=getarcclass(ocStruct,'lagrangemultcc');
            varStruct.(varname).arcidentifier=field2arcidentifier(arcid);
            for ii=1:numel(arcid)
                lagrangemultcc=algebraicterm2string(ocStruct,'symboliclagrangemultcc',0,arcid{ii});
                for jj=1:numel(lagrangemultcc.term)
                    if jj==1
                        varStruct.(varname).string{ii}{jj}=['out=' lagrangemultcc.term{jj}];
                    else
                        varStruct.(varname).string{ii}{jj}=lagrangemultcc.term{jj};
                    end
                end
                varStruct.(varname).string{ii}{jj}=[varStruct.(varname).string{ii}{jj} ';'];
            end
            varStruct.(varname).arcdependent=1;
            varStruct.(varname).multline=1;
            varStruct.(varname).vectorize=0;
            varStruct.(varname).type='algterm';
        catch
            %lasterr
            varStruct.(varname).string=varname;
        end

    case 'SYMBOLICSTATELAGRANGEMULTIPLIER'
        try
            arcid=getarcclass(ocStruct,'lagrangemultsc');
            varStruct.(varname).arcidentifier=field2arcidentifier(arcid);
            for ii=1:numel(arcid)
                lagrangemultsc=algebraicterm2string(ocStruct,'symboliclagrangemultsc',0,arcid{ii});
                for jj=1:numel(lagrangemultsc.term)
                    if jj==1
                        varStruct.(varname).string{ii}{jj}=['out=' lagrangemultsc.term{jj}];
                    else
                        varStruct.(varname).string{ii}{jj}=lagrangemultsc.term{jj};
                    end
                end
                varStruct.(varname).string{ii}{jj}=[varStruct.(varname).string{ii}{jj} ';'];
            end
            varStruct.(varname).arcdependent=1;
            varStruct.(varname).multline=1;
            varStruct.(varname).vectorize=0;
            varStruct.(varname).type='algterm';
        catch
            %lasterr
            varStruct.(varname).string=varname;
        end

    case 'ARCIDENTIFIERCONTROL'
        varStruct.(varname).arcdependent=1;
        try
            arcfield=getarcclass(ocStruct,'control');
            for ii=1:numel(arcfield)
                arcidentifier=field2arcidentifier(arcfield{ii});
                if numel(arcidentifier)>1
                    varStruct.(varname).string{ii}=['{' arcidentifier{1}];
                    for jj=2:numel(arcidentifier)
                        varStruct.(varname).string{ii}=[varStruct.(varname).string{ii} ',' arcidentifier{jj}];
                    end
                    varStruct.(varname).string{ii}=[varStruct.(varname).string{ii} '}'];
                else
                    varStruct.(varname).string{ii}=arcidentifier{1};
                end
            end
        catch
            %lasterr
            varStruct.(varname).string={'0'};
        end
    case 'SYMBOLICDISCIMPULSEOBJECTIVEFUNCTION'
        try
            objectivefunction=algebraicterm2string(ocStruct,'symbolicdiscimpulseobjectivefunction',0);
            varStruct.(varname).string{1}=['out=' objectivefunction.term{1} ';'];
            varStruct.(varname).multline=1;
            varStruct.(varname).vectorize=0;
            varStruct.(varname).type='algterm';
        catch
            %lasterr
            varStruct.(varname).string=varname;
        end

    case 'SYMBOLICDISCOBJECTIVEFUNCTION'
        try
            objectivefunction=algebraicterm2string(ocStruct,'symbolicdiscobjectivefunction',0);
            varStruct.(varname).string{1}=['out=' objectivefunction.term{1} ';'];
            varStruct.(varname).multline=1;
            varStruct.(varname).vectorize=0;
            varStruct.(varname).type='algterm';
        catch
            %lasterr
            varStruct.(varname).string=varname;
        end

    case 'DISCOBJECTIVEFUNCTION'
        try
            objectivefunction=algebraicterm2string(ocStruct,'discobjectivefunction',1);
            varStruct.(varname).string{1}=['out=' objectivefunction.term{1} ';'];
            varStruct.(varname).multline=1;
            varStruct.(varname).vectorize=1;
            varStruct.(varname).type='algterm';
        catch
            %lasterr
            varStruct.(varname).string=varname;
        end

    case 'DISCOBJECTIVEFUNCTIONDERIVATIVETIME'
        try
            arcfield=getarcclass(ocStruct,'canonicalsystemjacobian');
            varStruct.(varname).arcidentifier=field2arcidentifier(arcfield);
            for ii=1:numel(arcfield)
                discobjectivefunctionderivativetime=algebraicterm2string(ocStruct,'discobjectivefunctionderivativetime',0,arcfield{ii});
                for jj=1:numel(discobjectivefunctionderivativetime.term)
                    if jj==1
                        varStruct.(varname).string{ii}{jj}=['out=' discobjectivefunctionderivativetime.term{jj}];
                    else
                        varStruct.(varname).string{ii}{jj}=discobjectivefunctionderivativetime.term{jj};
                    end
                end
                varStruct.(varname).string{ii}{jj}=[varStruct.(varname).string{ii}{jj} ';'];
            end
            varStruct.(varname).multline=1;
            arcnum=retrieveimpulsemodelinformation(ocStruct,'arcnum');
            if arcnum.value>1
                varStruct.(varname).arcdependent=1;
            else
                varStruct.(varname).arcdependent=0;
                varStruct.(varname).string=varStruct.(varname).string{1};
            end
            varStruct.(varname).type='algterm';
        catch
            %lasterr
            varStruct.(varname).string=varname;
        end

    case 'DISCOBJECTIVEFUNCTIONJACOBIAN'
        try
            arcfield=getarcclass(ocStruct,'canonicalsystemjacobian');
            arcnum=retrieveimpulsemodelinformation(ocStruct,'arcnum');
            varStruct.(varname).arcidentifier=field2arcidentifier(arcfield);
            for ii=1:arcnum.value
                discobjectivefunctionDx=algebraicterm2string(ocStruct,'discobjectivefunctionjacobian',0,arcfield{ii});
                for jj=1:numel(discobjectivefunctionDx.term)
                    if jj==1
                        varStruct.(varname).string{ii}{jj}=['out=' discobjectivefunctionDx.term{jj}];
                    else
                        varStruct.(varname).string{ii}{jj}=discobjectivefunctionDx.term{jj};
                    end
                end
                varStruct.(varname).string{ii}{jj}=[varStruct.(varname).string{ii}{jj} ';'];
            end
            varStruct.(varname).multline=1;
            varStruct.(varname).vectorize=0;
            varStruct.(varname).type='algterm';
            if arcnum.value>1
                varStruct.(varname).arcdependent=1;
            else
                varStruct.(varname).arcdependent=0;
                varStruct.(varname).string=varStruct.(varname).string{1};
            end
        catch
            %lasterr
            varStruct.(varname).string=varname;
        end

    case 'DISCOBJECTIVEFUNCTIONPARAMETERJACOBIAN'
        try
            arcfield=getarcclass(ocStruct,'canonicalsystemjacobian');
            arcnum=retrieveimpulsemodelinformation(ocStruct,'arcnum');
            varStruct.(varname).arcidentifier=field2arcidentifier(arcfield);
            for ii=1:arcnum.value
                discobjectivefunctionDx=algebraicterm2string(ocStruct,'discobjectivefunctionparameterjacobian',0,arcfield{ii});
                for jj=1:numel(discobjectivefunctionDx.term)
                    if jj==1
                        varStruct.(varname).string{ii}{jj}=['out=' discobjectivefunctionDx.term{jj}];
                    else
                        varStruct.(varname).string{ii}{jj}=discobjectivefunctionDx.term{jj};
                    end
                end
                varStruct.(varname).string{ii}{jj}=[varStruct.(varname).string{ii}{jj} ';'];
            end
            varStruct.(varname).multline=1;
            varStruct.(varname).vectorize=0;
            varStruct.(varname).type='algterm';
            if arcnum.value>1
                varStruct.(varname).arcdependent=1;
            else
                varStruct.(varname).arcdependent=0;
                varStruct.(varname).string=varStruct.(varname).string{1};
            end
        catch
            %lasterr
            varStruct.(varname).string=varname;
        end

    case 'DISCIMPULSEOBJECTIVEFUNCTION'
        try
            objectivefunction=algebraicterm2string(ocStruct,'discimpulseobjectivefunction',1);
            varStruct.(varname).string{1}=['out=' objectivefunction.term{1} ';'];
            varStruct.(varname).multline=1;
            varStruct.(varname).vectorize=1;
            varStruct.(varname).type='algterm';
        catch
            %lasterr
            varStruct.(varname).string=varname;
        end
        
    case 'HAMILTONIANFUNCTION'
        try
            hamiltonianfunction=algebraicterm2string(ocStruct,'hamiltonianfunction',1);
            varStruct.(varname).string{1}=['out=' hamiltonianfunction.term{1} ';'];
            varStruct.(varname).multline=1;
            varStruct.(varname).vectorize=1;
            varStruct.(varname).type='algterm';
        catch
            %lasterr
            varStruct.(varname).string=varname;
        end

    case 'IMPULSEHAMILTONIANFUNCTION'
        try
            hamiltonianfunction=algebraicterm2string(ocStruct,'impulsehamiltonianfunction',1);
            varStruct.(varname).string{1}=['out=' hamiltonianfunction.term{1} ';'];
            varStruct.(varname).multline=1;
            varStruct.(varname).vectorize=0;
            varStruct.(varname).type='algterm';
        catch
            %lasterr
            varStruct.(varname).string=varname;
        end


    case 'SALVAGEVALUE'
        try
            salvagevalue=algebraicterm2string(ocStruct,'salvagevalue',1);
            varStruct.(varname).string{1}=['out=' salvagevalue.term ';'];
            varStruct.(varname).multline=0;
            varStruct.(varname).vectorize=0;
            varStruct.(varname).type='algterm';
        catch
            %lasterr
            varStruct.(varname).string=varname;
        end

    case 'SYMBOLICSALVAGEVALUE'
        try
            salvagevalue=algebraicterm2string(ocStruct,'symbolicsalvagevalue',0);
            varStruct.(varname).string{1}=['out=' salvagevalue.term{1} ';'];
            varStruct.(varname).multline=0;
            varStruct.(varname).vectorize=0;
            varStruct.(varname).type='algterm';
        catch
            %lasterr
            varStruct.(varname).string=varname;
        end

    case 'DISCOUNTEDSALVAGEVALUE'
        try
            salvagevalue=algebraicterm2string(ocStruct,'discountedsalvagevalue',1);
            varStruct.(varname).string{1}=['out=' salvagevalue.term ';'];
            varStruct.(varname).multline=0;
            varStruct.(varname).vectorize=0;
            varStruct.(varname).type='algterm';
        catch
            %lasterr
            varStruct.(varname).string=varname;
        end

    case 'IMPULSEOBJECTIVEFUNCTION'
        try
            objectivesummand=algebraicterm2string(ocStruct,'objectivesummand',1);
            varStruct.(varname).string{1}=['out=' objectivesummand.term ';'];
            varStruct.(varname).multline=0;
            varStruct.(varname).vectorize=0;
            varStruct.(varname).type='algterm';
        catch
            %lasterr
            varStruct.(varname).string=varname;
        end

    case 'DISCOUNTEDIMPULSEOBJECTIVEFUNCTION'
        try
            objectivesummand=algebraicterm2string(ocStruct,'discountedobjectivesummand',1);
            varStruct.(varname).string{1}=['out=' objectivesummand.term ';'];
            varStruct.(varname).multline=0;
            varStruct.(varname).vectorize=0;
            varStruct.(varname).type='algterm';
        catch
            %lasterr
            varStruct.(varname).string=varname;
        end
    case 'TRANSVERSALITYBC'
        try
            DsalvagevalueDx=algebraicterm2string(ocStruct,'DsalvagevalueDx',1);
            for ii=1:numel(DsalvagevalueDx.term)
                if ii==1
                    varStruct.(varname).string{ii}=['out=' DsalvagevalueDx.term{ii}];
                else
                    varStruct.(varname).string{ii}=DsalvagevalueDx.term{ii};
                end
            end
            varStruct.(varname).string{ii}=[varStruct.(varname).string{ii} ';'];
            varStruct.(varname).arcdependent=0;
            varStruct.(varname).multline=1;
            varStruct.(varname).vectorize=0;
            varStruct.(varname).type='algterm';
        catch
            %lasterr
            varStruct.(varname).string=varname;
        end


    case 'PONTRYAGINFUNCTION'
        try
            pontryaginfunction=algebraicterm2string(ocStruct,'pontryaginfunction',1);
            varStruct.(varname).string{1}=['out=' pontryaginfunction.term{1} ';'];
            varStruct.(varname).multline=1;
            varStruct.(varname).vectorize=1;
            varStruct.(varname).type='algterm';
        catch
            %lasterr
            varStruct.(varname).string=varname;
        end


    case 'SYMBOLICOBJECTIVEFUNCTION'
        try
            symbolicobjectivefunction=algebraicterm2string(ocStruct,'symbolicobjectivefunction',0);
            varStruct.(varname).string{1}=['out=' symbolicobjectivefunction.term{1} ';'];
            varStruct.(varname).multline=1;
            varStruct.(varname).vectorize=0;
            varStruct.(varname).type='algterm';
        catch
            %lasterr
            varStruct.(varname).string=varname;
        end

    case 'SYMBOLICPONTRYAGINFUNCTION'
        try
            symbolicpontryaginfunction=algebraicterm2string(ocStruct,'symbolicpontryaginfunction',0);
            varStruct.(varname).string{1}=['out=' symbolicpontryaginfunction.term{1} ';'];
            varStruct.(varname).multline=1;
            varStruct.(varname).vectorize=0;
            varStruct.(varname).type='algterm';
        catch
            %lasterr
            varStruct.(varname).string=varname;
        end

    case 'SYMBOLICIMPULSEPONTRYAGINFUNCTION'
        try
            symbolicpontryaginfunction=algebraicterm2string(ocStruct,'symbolicimpulsepontryaginfunction',0);
            varStruct.(varname).string{1}=['out=' symbolicpontryaginfunction.term{1} ';'];
            varStruct.(varname).multline=1;
            varStruct.(varname).vectorize=0;
            varStruct.(varname).type='algterm';
        catch
            %lasterr
            varStruct.(varname).string=varname;
        end

    case 'DHAMILTONIANDX'
        try
            arcfield=getarcclass(ocStruct,'canonicalsystemjacobian');
            arcnum=retrieveimpulsemodelinformation(ocStruct,'arcnum');
            varStruct.(varname).arcidentifier=field2arcidentifier(arcfield);
            for ii=1:arcnum.value
                pontryaginfunctionDx=algebraicterm2string(ocStruct,'pontryaginfunctionDX',0,arcfield{ii});
                for jj=1:numel(pontryaginfunctionDx.term)
                    if jj==1
                        varStruct.(varname).string{ii}{jj}=['out=' pontryaginfunctionDx.term{jj}];
                    else
                        varStruct.(varname).string{ii}{jj}=pontryaginfunctionDx.term{jj};
                    end
                end
                varStruct.(varname).string{ii}{jj}=[varStruct.(varname).string{ii}{jj} ';'];
            end
            varStruct.(varname).multline=1;
            varStruct.(varname).vectorize=0;
            varStruct.(varname).type='algterm';
            if arcnum.value>1
                varStruct.(varname).arcdependent=1;
            else
                varStruct.(varname).arcdependent=0;
                varStruct.(varname).string=varStruct.(varname).string{1};
            end
        catch
            %lasterr
            varStruct.(varname).string=varname;
        end
    case 'DIMPULSEHAMILTONIANDTAU'
        try
            impulsepontryaginfunctionDtau=algebraicterm2string(ocStruct,'impulsepontryaginfunctionDtau',0);
                    varStruct.(varname).string{1}=['out=' impulsepontryaginfunctionDtau.term{1} ';'];
            varStruct.(varname).multline=1;
            varStruct.(varname).vectorize=0;
            varStruct.(varname).type='algterm';
            varStruct.(varname).arcdependent=0;
        catch
            %lasterr
            varStruct.(varname).string=varname;
        end

    case 'DATAPATHNAME'
        varStruct.(varname).string=['''' fullocmatfile(getocmatfolder('userdata',modeltype(ocStruct),modelname(ocStruct))) ''''];

    case 'RESULTFILE'
        varStruct.(varname).string=['''' fullocmatfile(getocmatfolder('userdata',modeltype(ocStruct),modelname(ocStruct)),'SaveIntermediateResults') ''''];

    case 'GLOBALVARFILE'
        varStruct.(varname).string=['''' fullocmatfile(getocmatfolder('userdata',modeltype(ocStruct),modelname(ocStruct)),'SaveIntermediateResultsGlobalVariable') ''''];

    case 'INFODETAILS'
        varStruct.(varname).multline=2;
        varStruct.(varname).string={'%', ...
            ['% this file was automatically created: ' datestr(now)], ...
            '% written by Dieter Grass, 2013'};
        
%% for the gradient method        
    case 'GRADIENTHAMILTONIAN'
        try
            hamiltonianDu=algebraicterm2string(ocStruct,'hamiltonianDu',1);
            for jj=1:numel(hamiltonianDu.term)
                if jj==1
                    varStruct.(varname).string{jj}=['out=' hamiltonianDu.term{jj}];
                else
                    varStruct.(varname).string{jj}=hamiltonianDu.term{jj};
                end
            end
            varStruct.(varname).string{jj}=[varStruct.(varname).string{jj} ';'];
            varStruct.(varname).multline=1;
            varStruct.(varname).vectorize=1;
            varStruct.(varname).type='algterm';
            varStruct.(varname).arcdependent=0;
        catch
            %lasterr
            varStruct.(varname).string=varname;
        end

    case 'STATEDYNAMICS'
        try
            statedynamics=algebraicterm2string(ocStruct,'statedynamics',1,'');
            for jj=1:numel(statedynamics.term)
                if jj==1
                    varStruct.(varname).string{jj}=['out=' statedynamics.term{jj}];
                else
                    varStruct.(varname).string{jj}=statedynamics.term{jj};
                end
            end
            varStruct.(varname).string{jj}=[varStruct.(varname).string{jj} ';'];
            %end
            varStruct.(varname).multline=1;
            varStruct.(varname).vectorize=1;
            varStruct.(varname).arcdependent=0;
            varStruct.(varname).type='algterm';
        catch
            %lasterr
            varStruct.(varname).string=varname;
        end

    case 'ADJOINTDYNAMICS'
        try
            adjointsystem=algebraicterm2string(ocStruct,'adjointsystem',1,'');
            for jj=1:numel(adjointsystem.term)
                if jj==1
                    varStruct.(varname).string{jj}=['out=' adjointsystem.term{jj}];
                else
                    varStruct.(varname).string{jj}=adjointsystem.term{jj};
                end
            end
            varStruct.(varname).string{jj}=[varStruct.(varname).string{jj} ';'];
            %end
            varStruct.(varname).multline=1;
            varStruct.(varname).vectorize=1;
            varStruct.(varname).arcdependent=0;
            varStruct.(varname).type='algterm';
        catch
            %lasterr
            varStruct.(varname).string=varname;
        end
    case 'STATEEVENTVALUE'
        try
            stateevent=algebraicterm2string(ocStruct,'stateevent',1,'');
            for jj=1:numel(stateevent.term)
                if jj==1
                    varStruct.(varname).string{jj}=['out=' stateevent.term{jj}];
                else
                    varStruct.(varname).string{jj}=stateevent.term{jj};
                end
            end
            varStruct.(varname).string{jj}=[varStruct.(varname).string{jj} ';'];
            %end
            varStruct.(varname).multline=1;
            varStruct.(varname).vectorize=0;
            varStruct.(varname).arcdependent=0;
            varStruct.(varname).type='algterm';
        catch
            %lasterr
            varStruct.(varname).string=varname;
        end

    case 'ADJOINTEVENTVALUE'
        try
            adjointevent=algebraicterm2string(ocStruct,'adjointevent',1,'');
            for jj=1:numel(adjointevent.term)
                if jj==1
                    varStruct.(varname).string{jj}=['out=' adjointevent.term{jj}];
                else
                    varStruct.(varname).string{jj}=adjointevent.term{jj};
                end
            end
            varStruct.(varname).string{jj}=[varStruct.(varname).string{jj} ';'];
            %end
            varStruct.(varname).multline=1;
            varStruct.(varname).vectorize=0;
            varStruct.(varname).arcdependent=0;
            varStruct.(varname).type='algterm';
        catch
            %lasterr
            varStruct.(varname).string=varname;
        end
        
    case 'OBJECTIVEFUNCTION'
        try
            objectivefunction=algebraicterm2string(ocStruct,'objectivefunction',1);
            varStruct.(varname).string{1}=['out=' objectivefunction.term{1} ';'];
            varStruct.(varname).multline=0;
            varStruct.(varname).vectorize=1;
            varStruct.(varname).type='algterm';
        catch
            %lasterr
            varStruct.(varname).string=varname;
        end


    case 'TRANSVERSALITYCONDITION'
        try
            transversalitycondition=algebraicterm2string(ocStruct,'transversalitycondition',0);
            for jj=1:numel(transversalitycondition.term)
                if jj==1
                    varStruct.(varname).string{jj}=['out=' transversalitycondition.term{jj}];
                else
                    varStruct.(varname).string{jj}=transversalitycondition.term{jj};
                end
            end
            varStruct.(varname).string{jj}=[varStruct.(varname).string{jj} ';'];
            varStruct.(varname).multline=1;
            varStruct.(varname).vectorize=0;
            varStruct.(varname).type='algterm';
        catch
            %lasterr
            varStruct.(varname).string=varname;
        end

end

function out=algebraicterm2string(varargin)

out=impulsemodelalgebraicterm2string(varargin{:});