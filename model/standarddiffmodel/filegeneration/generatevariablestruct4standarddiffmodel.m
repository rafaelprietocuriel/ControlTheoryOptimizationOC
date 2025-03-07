function varStruct=generatevariablestruct4standarddiffmodel(ocStruct,varStruct,varname,overwrite,varargin)
%
% a structure is generated which is then used to replace variables in
% template model files
% with overwrite 1 an already existing field of varStruct is overwritten
% the structure consists of a field given by the name of the variable. If
% the actual value of the variable depends on the arc the field string
% contains a cell array of strings and the field arcependent is set to one.
% multline 1 ... lines>1 are intended, e.g., DYNAMICSDYNAMICS
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

try
    switch varname
        case 'INDEPENDENT'
            varStruct.(varname).string=getbasicname('independent');

        case 'DEPENDENTVAR'
            varStruct.(varname).string=getbasicname('dependent');

        case 'CONTROL'
            varStruct.(varname).string=getbasicname('control');

        case 'PARVAR'
            varStruct.(varname).string=getbasicname('parametervariables');

        case 'ARCVAR'
            varStruct.(varname).string=getbasicname('arcidentifiervar');

        case 'ARCIDENTIFIER'
            varStruct.(varname).arcdependent=1;
            arcidentifier=retrievediffmodelinformation(ocStruct,'arcidentifier');
            varStruct.(varname).string=arcidentifier.value;

        case 'ODEDIM'
            varStruct.(varname).arcdependent=1;
            odedim=retrievediffmodelinformation(ocStruct,'odedim');
            arcnum=retrievediffmodelinformation(ocStruct,'arcnum');
            for ii=1:arcnum.value
                varStruct.(varname).string{ii}=int2str(odedim.value);
            end

        case 'AEDIM'
            varStruct.(varname).arcdependent=1;
            arcnum=retrievediffmodelinformation(ocStruct,'arcnum');
            arcidentifier=retrievediffmodelinformation(ocStruct,'arcidentifier');
            for ii=1:arcnum.value
                algebraicequationnum=retrievediffmodelinformation(ocStruct,'algebraicequationnum',arcidentifier.value{ii});
                varStruct.(varname).string{ii}=int2str(algebraicequationnum.value);
            end

        case 'DAEORDER'
            if ~isfield(varStruct,'AEDIM')
                varStruct=generatevariablestruct4standardmodel(ocStruct,varStruct,'AEDIM');
            end
            varStruct.(varname).arcdependent=1;
            arcnum=retrievediffmodelinformation(ocStruct,'arcnum');
            statenum=retrievediffmodelinformation(ocStruct,'statenum');
            for ii=1:arcnum.value
                order=['[' int2str(ones(1,statenum.value))];
                aedim=str2double(varStruct.AEDIM.string{ii});
                if ~isnan(aedim) && aedim
                    order=[order '  ' int2str(zeros(1, aedim))];
                end
                order=[order ']'];
                varStruct.(varname).string{ii}=order;
            end

        case 'SUBSFLAG'
            varStruct.(varname).string='s';

        case 'EXOGENOUSFUNCTION'
            varStruct.(varname).string=getbasicname('exogenousfunction');


        case 'CONSTRAINTNUM'
            inequalitycontrolconstraintnum=retrievediffmodelinformation(ocStruct,'inequalitycontrolconstraintnum');
            inequalitystateconstraintnum=retrievediffmodelinformation(ocStruct,'inequalitystateconstraintnum');
            varStruct.(varname).string=num2str(inequalitycontrolconstraintnum.value+inequalitystateconstraintnum.value);


        case 'CONTROLCONSTRAINTNUM'
            inequalitycontrolconstraintnum=retrievediffmodelinformation(ocStruct,'inequalitycontrolconstraintnum');
            varStruct.(varname).string='';
            varStruct.(varname).string=num2str(inequalitycontrolconstraintnum.value);

        case 'STATECONSTRAINTNUM'
            inequalitystateconstraintnum=retrievediffmodelinformation(ocStruct,'inequalitystateconstraintnum');
            varStruct.(varname).string='';
            varStruct.(varname).string=num2str(inequalitystateconstraintnum.value);

        case 'COSTATECOORD'
            statenum=retrievediffmodelinformation(ocStruct,'statenum');
            costatenum=retrievediffmodelinformation(ocStruct,'costatenum');
            varStruct.(varname).string=[int2str(statenum.value+1) ':' int2str(statenum.value+costatenum.value)];

        case 'ODIFFDIM'
            odedim=retrievediffmodelinformation(ocStruct,'odedim');
            varStruct.(varname).string=int2str(odedim.value);

        case 'ODIFFDIMCOORD'
            odedim=retrievediffmodelinformation(ocStruct,'odedim');
            varStruct.(varname).string=['1:' int2str(odedim.value)];

        case 'LAGRANGEMULTCC'
            varStruct.(varname).string=getbasicname('lagrangemultcc');

        case 'LAGRANGEMULTSC'
            varStruct.(varname).string=getbasicname('lagrangemultsc');

        case 'GENERATIONDATE'
            varStruct.(varname).string=datestr(now);

        case 'MODELNAME'
            varStruct.(varname).string=modelname(ocStruct);

        case 'MODELTYPE'
            varStruct.(varname).string=modeltype(ocStruct);

        case 'NUMBEROFARCS'
            arcnum=retrievediffmodelinformation(ocStruct,'arcnum');
            varStruct.(varname).string=num2str(arcnum.value);

        case 'ARCCOLOR'
            arcnum=retrievediffmodelinformation(ocStruct,'arcnum');
            defaultcolor=get(0,'DefaultAxesColorOrder');
            defaultcolor=repmat(defaultcolor,ceil(arcnum.value/size(defaultcolor,1)),1);
            varStruct.(varname).string=mat2str(defaultcolor(1:arcnum.value,:));
        case 'ARCARGUMENT'
            arcargument=retrievediffmodelinformation(ocStruct,'argument');
            varStruct.(varname).string=mat2str(arcargument.value);

        case 'EDGES'
            arcargument=retrievediffmodelinformation(ocStruct,'argument');
            arcnum=retrievediffmodelinformation(ocStruct,'arcnum');
            if arcnum.value>1
                nc=nchoosek(arcargument.value,2);
                varStruct.(varname).string=mat2str([nc;nc(:,[2 1])]);
            else
                varStruct.(varname).string='[]';
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

        case 'ZEROSNTIMESNUMEQ'
            % returns a zero matrix where number of rows = number of states and
            % number of columns = number of equations of the canonical system
            try
                arcfield=getarcclass(ocStruct,'algebraicequation');
                statenum=retrievediffmodelinformation(ocStruct,'statenum');
                for ii=1:numel(arcfield)
                    arcidentifier=field2arcidentifier(arcfield{ii});
                    equationnum=retrievediffmodelinformation(ocStruct,'canonicalsystemequationnum',arcidentifier{1});
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

        case 'ISAUTONOMOUS'
            autonomousflag=retrievediffmodelinformation(ocStruct,'autonomous');
            varStruct.(varname).string=num2str(autonomousflag.value);

        case 'EXOGENOUSFUNCTIONNUM'
            exogenousfunctionnum=retrievediffmodelinformation(ocStruct,'exogenousfunctionnum');
            varStruct.(varname).string=num2str(exogenousfunctionnum.value);

        case 'STATECOORD'
            statenum=retrievediffmodelinformation(ocStruct,'statenum');
            varStruct.(varname).string=['1:' int2str(statenum.value)];

        case 'STATENUM'
            statenum=retrievediffmodelinformation(ocStruct,'statenum');
            varStruct.(varname).string=int2str(statenum.value);

        case 'ZEROSONESTATENUM'
            statenum=retrievediffmodelinformation(ocStruct,'statenum');
            varStruct.(varname).string=mat2str(zeros(1,statenum.value));

        case 'ZEROSTWOSTATENUM'
            statenum=retrievediffmodelinformation(ocStruct,'statenum');
            varStruct.(varname).string=mat2str(zeros(2,statenum.value));


        case 'STATECOSTATECOORD'
            statenum=retrievediffmodelinformation(ocStruct,'statenum');
            costatenum=retrievediffmodelinformation(ocStruct,'costatenum');
            varStruct.(varname).string=['1:' int2str(statenum.value+costatenum.value)];
        case 'DEFAULTPLOTCOMMAND'
            depvarname=getbasicname('dependent');
            varStruct.(varname).string=['h=plot(' depvarname '(1,:),' depvarname '(2,:),''x'');'];

        case 'DEFAULTPLOTTHREECOMMAND'
            depvarname=getbasicname('dependent');
            varStruct.(varname).string=['h=plot3(' depvarname '(1,:),' depvarname '(2,:),' depvarname '(3,:),''x'');'];

        case 'IMPLICITCONTROLS'
            arcidentifier=retrievediffmodelinformation(ocStruct,'arcidentifier');
            b=false;
            for ii=1:length(arcidentifier.value)
                out=retrievediffmodelinformation(ocStruct,'implicitnonlinearcontrolindex',arcidentifier.value{ii});
                if ~isempty(out.value)
                    b=true;
                    return
                end
            end
            varStruct.(varname).string=num2str(b);
        case 'PARAMETERVALUES'
            % the strings used for assigning parametervalues from the par
            % argument in functions to the spcific paraemter variable names,
            % e.g. r=pararg(1); ...
            varStruct.(varname).multline=2;
            basename=getbasicname('parametervariables');
            parvar=algebraicodeterm2string(ocStruct,'parametervariables');
            for ii=1:numel(parvar.term)
                varStruct.(varname).string{ii}=[parvar.term{ii} '=' basename '(' num2str(ii) ');'];
            end

        case 'LATEXUSERDEPENDENTNAMEFIRST'
            statename=retrievediffmodelinformation(ocStruct,'statename');
            [val1 val2 val3]=regexp(statename.value{1},'([0-9]+)\>','split');
            if ~isempty(val2)
                varStruct.(varname).string=['$' latex(sym(val1{1})) '_' statename.value{1}(val2:val3) '$'];
            else
                varStruct.(varname).string=['$' latex(sym(val1{1})) '$'];
            end

        case 'LATEXUSERDEPENDENTNAMESECOND'
            statenum=retrievediffmodelinformation(ocStruct,'statenum');
            if statenum.value>1
                statename=retrievediffmodelinformation(ocStruct,'statename');
                labelname=statename.value{2};
            else
                costatename=retrievediffmodelinformation(ocStruct,'costatename');
                labelname=costatename.value{1};
            end
            [val1 val2 val3]=regexp(labelname,'([0-9]+)\>','split');
            if ~isempty(val2)
                varStruct.(varname).string=['$' latex(sym(val1{1})) '_' labelname(val2:val3) '$'];
            else
                varStruct.(varname).string=['$' latex(sym(val1{1})) '$'];
            end

        case 'LATEXUSERDEPENDENTNAMETHIRD'
            statenum=retrievediffmodelinformation(ocStruct,'statenum');
            if statenum.value<=2
                labelname='';
            else
                statename=retrievediffmodelinformation(ocStruct,'statename');
                labelname=statename.value{3};
            end

            [val1 val2 val3]=regexp(labelname,'([0-9]+)\>','split');
            if ~isempty(val2)
                varStruct.(varname).string=['$' latex(sym(val1{1})) '_' labelname(val2:val3) '$'];
            else
                varStruct.(varname).string=['$' latex(sym(val1{1})) '$'];
            end

        case 'EXOGENOUSFUNCTIONTERM'
            exogenousfunctionterm=algebraicodeterm2string(ocStruct,'exogenousfunctionterm',1,'');
            for jj=1:numel(exogenousfunctionterm.term)
                if jj==1
                    varStruct.(varname).string{jj}=exogenousfunctionterm.term{jj};
                else
                    varStruct.(varname).string{jj}=exogenousfunctionterm.term{jj};
                end
            end
            varStruct.(varname).string{jj}=[varStruct.(varname).string{jj} ';'];
            varStruct.(varname).multline=1;
            varStruct.(varname).vectorize=1;
            varStruct.(varname).arcdependent=0;
            varStruct.(varname).type='algterm';

        case 'EXPLICITDIFFERENCEEQUATION'
            explicitdifferenceequationflag=retrievediffmodelinformation(ocStruct,'maptype');
            varStruct.(varname).string=num2str(strcmp(explicitdifferenceequationflag.value,'explicit'));

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
        case 'CONTROLLAGRANGEMULTIPLIER'
            try
                arcid=getarcclass(ocStruct,'lagrangemultcc');
                varStruct.(varname).arcidentifier=field2arcidentifier(arcid);
                for ii=1:numel(arcid)
                    lagrangemultcc=algebraicstandarddiffmodelterm2string(ocStruct,'lagrangemultcc',1,arcid{ii});
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

        case 'CONSTRAINT'
            try
                constraint=algebraicstandarddiffmodelterm2string(ocStruct,'constraint',1);
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
        case 'STATELAGRANGEMULTIPLIER'
            try
                arcid=getarcclass(ocStruct,'lagrangemultsc');
                varStruct.(varname).arcidentifier=field2arcidentifier(arcid);
                if isempty(arcid)
                    varStruct.(varname).string='';
                else
                    for ii=1:numel(arcid)
                        lagrangemultsc=algebraicstandarddiffmodelterm2string(ocStruct,'lagrangemultsc',1,arcid{ii});
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
                    lagrangemultcc=algebraicstandarddiffmodelterm2string(ocStruct,'symboliclagrangemultcc',0,arcid{ii});
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
                    lagrangemultsc=algebraicstandarddiffmodelterm2string(ocStruct,'symboliclagrangemultsc',1,arcid{ii});
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
        case 'ADJOINTSTATETPONE'
            try
                arcfield=getarcclass(ocStruct,'canonicalsystemjacobian');
                arcnum=retrievediffmodelinformation(ocStruct,'arcnum');
                varStruct.(varname).arcidentifier=field2arcidentifier(arcfield);
                for ii=1:numel(arcfield)
                    adjointstateexplicit=algebraicstandarddiffmodelterm2string(ocStruct,'adjointstateexplicit',1,arcfield{ii});
                    for jj=1:numel(adjointstateexplicit.term)
                        if jj==1
                            varStruct.(varname).string{ii}{jj}=['costatetp1=' adjointstateexplicit.term{jj}];
                        else
                            varStruct.(varname).string{ii}{jj}=adjointstateexplicit.term{jj};
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

        case 'CANONICALSYSTEMMAP'
            try
                canonicalsystem=algebraicstandarddiffmodelterm2string(ocStruct,'canonicalsystem',1,'');
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


        case 'SYMBOLICCANONICALSYSTEMMAP'
            try
                canonicalsystem=algebraicstandarddiffmodelterm2string(ocStruct,'symboliccanonicalsystemmap',0,'');
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

        case 'CONTROLVALUESUBSTITUTION'
            controlnum=retrievediffmodelinformation(ocStruct,'controlnum');
            controlnametp1=retrievediffmodelinformation(ocStruct,'controlnametp1');
            inequalitycontrolconstraintnum=retrievediffmodelinformation(ocStruct,'inequalitycontrolconstraintnum');
            lagrangemultipliercontrolnametp1=retrievediffmodelinformation(ocStruct,'lagrangemultipliercontrolnametp1');
            indexstr='';
            for ii=1:controlnum.value
                indexstr=[indexstr '''' controlnametp1.value{ii} ''','];
            end
            coordstr=basename2vectorstring(getbasicname('control'),1:controlnum.value,'coord');
            if inequalitycontrolconstraintnum.value
                for ii=1:inequalitycontrolconstraintnum.value
                    indexstr=[indexstr '''' lagrangemultipliercontrolnametp1.value{ii} ''','];
                end
                %indexstr=[indexstr ''',''' basename2vectorstring(getbasicname('lagrangemultcc'),1:inequalitycontrolconstraintnum.value,'index',''',''')];
                coordstr=[coordstr ',' basename2vectorstring(getbasicname('lagrangemultcc'),1:inequalitycontrolconstraintnum.value,'coord')];
            end
            indexstr(end)=[];
            varStruct.(varname).string=['out=subs(out,{' indexstr '},{' coordstr '});'];
        case 'CANONICALSYSTEMMAPJACOBIAN'
            try
                arcfield=getarcclass(ocStruct,'canonicalsystemjacobian');
                arcnum=retrievediffmodelinformation(ocStruct,'arcnum');
                varStruct.(varname).arcidentifier=field2arcidentifier(arcfield);
                for ii=1:numel(arcfield)
                    canonicalsystemjacobian=algebraicstandarddiffmodelterm2string(ocStruct,'canonicalsystemjacobian',0,arcfield{ii});
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
                varStruct.(varname).vectorize=0;
            catch
                %lasterr
                varStruct.(varname).string=varname;
            end


        case 'SYMBOLICPONTRYAGINFUNCTION'
            try
                symbolicpontryaginfunction=algebraicstandarddiffmodelterm2string(ocStruct,'symbolicpontryaginfunction',0);
                varStruct.(varname).string{1}=['out=' symbolicpontryaginfunction.term{1} ';'];
                varStruct.(varname).multline=1;
                varStruct.(varname).vectorize=0;
                varStruct.(varname).type='algterm';
            catch
                %lasterr
                varStruct.(varname).string=varname;
            end


        case 'PONTRYAGINFUNCTION'
            try
                pontryaginfunction=algebraicstandarddiffmodelterm2string(ocStruct,'pontryaginfunction',1);
                varStruct.(varname).string{1}=['out=' pontryaginfunction.term{1} ';'];
                varStruct.(varname).multline=1;
                varStruct.(varname).vectorize=1;
                varStruct.(varname).type='algterm';
            catch
                %lasterr
                varStruct.(varname).string=varname;
            end

        case 'SYMBOLICCANONICALSYSTEMMAPJACOBIAN'
            try
                arcfield=getarcclass(ocStruct,'canonicalsystemjacobian');
                arcnum=retrievediffmodelinformation(ocStruct,'arcnum');
                varStruct.(varname).arcidentifier=field2arcidentifier(arcfield);
                for ii=1:numel(arcfield)
                    canonicalsystemjacobian=algebraicstandarddiffmodelterm2string(ocStruct,'symboliccanonicalsystemmapjacobian',0,arcfield{ii});
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

        case 'DISCOBJECTIVEFUNCTION'
            try
                objectivefunction=algebraicstandarddiffmodelterm2string(ocStruct,'discobjectivefunction',1);
                varStruct.(varname).string{1}=['out=' objectivefunction.term{1} ';'];
                varStruct.(varname).multline=1;
                varStruct.(varname).vectorize=1;
                varStruct.(varname).type='algterm';
            catch
                %lasterr
                varStruct.(varname).string=varname;
            end

        case 'CANONICALSYSTEMMAPPARAMETERJACOBIAN'
            try
                arcfield=getarcclass(ocStruct,'canonicalsystemjacobian');
                varStruct.(varname).arcidentifier=field2arcidentifier(arcfield);
                for ii=1:numel(arcfield)
                    canonicalsystemparameterjacobian=algebraicstandarddiffmodelterm2string(ocStruct,'canonicalsystemparameterjacobian',0,arcfield{ii});
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
                arcnum=retrievediffmodelinformation(ocStruct,'arcnum');
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

        case 'DISCOBJECTIVEFUNCTIONJACOBIANX'
            try
                arcfield=getarcclass(ocStruct,'canonicalsystemjacobian');
                arcnum=retrievediffmodelinformation(ocStruct,'arcnum');
                varStruct.(varname).arcidentifier=field2arcidentifier(arcfield);
                for ii=1:arcnum.value
                    discobjectivefunctionDx=algebraicstandarddiffmodelterm2string(ocStruct,'discobjectivefunctionjacobianx',1,arcfield{ii});
                    for jj=1:numel(discobjectivefunctionDx.term)
                        if jj==1
                            varStruct.(varname).string{ii}{jj}=['Jx=' discobjectivefunctionDx.term{jj}];
                        else
                            varStruct.(varname).string{ii}{jj}=discobjectivefunctionDx.term{jj};
                        end
                    end
                    varStruct.(varname).string{ii}{jj}=[varStruct.(varname).string{ii}{jj} ';'];
                end
                varStruct.(varname).multline=1;
                varStruct.(varname).vectorize=1;
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

        case 'DISCOBJECTIVEFUNCTIONJACOBIANL'
            try
                arcfield=getarcclass(ocStruct,'canonicalsystemjacobian');
                arcnum=retrievediffmodelinformation(ocStruct,'arcnum');
                varStruct.(varname).arcidentifier=field2arcidentifier(arcfield);
                for ii=1:arcnum.value
                    discobjectivefunctionDx=algebraicstandarddiffmodelterm2string(ocStruct,'discobjectivefunctionjacobianl',1,arcfield{ii});
                    for jj=1:numel(discobjectivefunctionDx.term)
                        if jj==1
                            varStruct.(varname).string{ii}{jj}=['Jl=' discobjectivefunctionDx.term{jj}];
                        else
                            varStruct.(varname).string{ii}{jj}=discobjectivefunctionDx.term{jj};
                        end
                    end
                    varStruct.(varname).string{ii}{jj}=[varStruct.(varname).string{ii}{jj} ';'];
                end
                varStruct.(varname).multline=1;
                varStruct.(varname).vectorize=1;
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
                arcnum=retrievediffmodelinformation(ocStruct,'arcnum');
                varStruct.(varname).arcidentifier=field2arcidentifier(arcfield);
                for ii=1:arcnum.value
                    discobjectivefunctionDx=algebraicstandarddiffmodelterm2string(ocStruct,'discobjectivefunctionparameterjacobian',0,arcfield{ii});
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

        case 'STATEDYNAMICS'
            try
                statedynamics=algebraicstandarddiffmodelterm2string(ocStruct,'statedynamics',1,'');
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

        case 'STATEDYNAMICSJACOBIAN'
            try
                arcfield=getarcclass(ocStruct,'canonicalsystemjacobian');
                arcnum=retrievediffmodelinformation(ocStruct,'arcnum');
                varStruct.(varname).arcidentifier=field2arcidentifier(arcfield);
                for ii=1:numel(arcfield)
                    specificstatedynamics=algebraicstandarddiffmodelterm2string(ocStruct,'specificstatedynamicsjacobian',0,arcfield{ii});
                    for jj=1:numel(specificstatedynamics.term)
                        if jj==1
                            varStruct.(varname).string{ii}{jj}=['out=' specificstatedynamics.term{jj}];
                        else
                            varStruct.(varname).string{ii}{jj}=specificstatedynamics.term{jj};
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
                varStruct.(varname).vectorize=0;
            catch
                %lasterr
                varStruct.(varname).string=varname;
            end

        case 'OPTIMALCONTROLVALUE'
            try
                arcfield=getarcclass(ocStruct,'controlvalue');
                varStruct.(varname).arcidentifier=field2arcidentifier(arcfield);
                arcnum=retrievediffmodelinformation(ocStruct,'arcnum');
                for ii=1:numel(arcfield)
                    controlvalue=algebraicstandarddiffmodelterm2string(ocStruct,'controlvalue',1,arcfield{ii});
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


        case 'SYMBOLICOPTIMALCONTROLVALUE'
            try
                arcfield=getarcclass(ocStruct,'controlvalue');
                varStruct.(varname).arcidentifier=field2arcidentifier(arcfield);
                arcnum=retrievediffmodelinformation(ocStruct,'arcnum');
                for ii=1:numel(arcfield)
                    controlvalue=algebraicstandarddiffmodelterm2string(ocStruct,'symboliccontrolvalue',0,arcfield{ii});
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
                if arcnum.value>1
                    varStruct.(varname).arcdependent=1;
                else
                    varStruct.(varname).arcdependent=0;
                    varStruct.(varname).string=varStruct.(varname).string{1};
                end
                varStruct.(varname).vectorize=0;
                varStruct.(varname).type='algterm';
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
    end
catch
    ocmatmsg('Problems with variable ''%s''.',varname)
    varStruct.(varname).string=varname;
end

