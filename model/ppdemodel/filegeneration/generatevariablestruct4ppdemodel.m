function varStruct=generatevariablestruct4ppdemodel(ocStruct,varStruct,varname,overwrite,varargin)
%
% a structure is generated which is then used to replace variables in
% template model files
% with overwrite 1 an already existing field of varStruct is overwritten
% the structure consists of a field given by the name of the variable. If
% the actual value of the variable depends on the arc the field string
% contains a cell array of strings and the field arcependent is set to one.
% multline 1 ... lines>1 are intended, e.g., CANONICALSYSTEMDYNAMICS
% multline 2 ...lines are not intended, e.g., PARAMETERVALUES
opt=[];
if nargin<=3
    overwrite=0;
end
if ~overwrite
    if isfield(varStruct,varname)
        return
    end
end
if nargin>=5
    opt=varargin{1};
end

varStruct.(varname).arcdependent=0;
varStruct.(varname).arcidentifier=[];
varStruct.(varname).vectorize=0;
varStruct.(varname).type='';
varStruct.(varname).multline=0;

switch varname
    % option dependent variables
    % model dependent variables
    case 'TIME'
        time=retrieveppdemodelinformation(ocStruct,'time');
        varStruct.(varname).string=time.value;
        
    case 'SPACE'
        space=retrieveppdemodelinformation(ocStruct,'space');
        varStruct.(varname).string=space.value;
        
    case 'SPACEMID'
        space=retrieveppdemodelinformation(ocStruct,'spacemid');
        varStruct.(varname).string=space.value;
        
    case 'MYSPACEARG'
        spacearg=retrieveppdemodelinformation(ocStruct,'space');
        varStruct.(varname).string=spacearg.value;
        
    case 'DEPENDENTVAR'
        varStruct.(varname).string=getbasicname('dependent');
        
        
    case 'JACOBIANEQUATIONNUM'
        try
            statenum=retrieveppdemodelinformation(ocStruct,'statenum');
            costatenum=retrieveppdemodelinformation(ocStruct,'costatenum');
            gridnum=getbasicname('femdatagridnum');
            varStruct.(varname).string=[int2str(statenum.value+costatenum.value) '*' gridnum];
            varStruct.(varname).arcdependent=0;
        catch
            %lasterr
            varStruct.(varname).string=varname;
        end
        
    case 'DEPENDENTVARMID'
        depvar=getbasicname('dependent');
        femdataXcoordm1=retrieveppdemodelinformation(ocStruct,'femdataXcoordm1');
        femdataXcoordp1=retrieveppdemodelinformation(ocStruct,'femdataXcoordp1');
        Xcoordm1='';
        Xcoordp1='';
        for ii=1:length(femdataXcoordm1.value)
            Xcoordm1=[Xcoordm1 ',' femdataXcoordm1.value{ii}];
            Xcoordp1=[Xcoordp1 ',' femdataXcoordp1.value{ii}];
        end
        Xcoordm1(1)=[];
        Xcoordp1(1)=[];
        varStruct.(varname).string=['(' depvar '([' Xcoordm1 '])+' depvar '([' Xcoordp1 ']))*0.5;'];
        
    case 'DEPENDENTVARMIDV'
        depvar=getbasicname('dependent');
        femdataXcoordm1=retrieveppdemodelinformation(ocStruct,'femdataXcoordm1');
        femdataXcoordp1=retrieveppdemodelinformation(ocStruct,'femdataXcoordp1');
        Xcoordm1='';
        Xcoordp1='';
        for ii=1:length(femdataXcoordm1.value)
            Xcoordm1=[Xcoordm1 ',' femdataXcoordm1.value{ii}];
            Xcoordp1=[Xcoordp1 ',' femdataXcoordp1.value{ii}];
        end
        Xcoordm1(1)=[];
        Xcoordp1(1)=[];
        varStruct.(varname).string=['(' depvar '([' Xcoordm1 '],:)+' depvar '([' Xcoordp1 '],:))*0.5;'];
        
    case 'ENDTIME'
        varStruct.(varname).string=getbasicname('endtime');
        
    case 'PARVAR'
        varStruct.(varname).string=getbasicname('parametervariables');
        
    case 'FEMDATA'
        varStruct.(varname).string=getbasicname('femdata');
        
    case 'FEMDATAGRID'
        femdatagrid=getbasicname('femdatagrid');
        femdatagrid=[femdatagrid '.x.'''];
        varStruct.(varname).string=femdatagrid;
        
    case 'FEMDATAGRIDMID'
        femdatagrid=getbasicname('femdatagridmid');
        femdatagrid=[femdatagrid '.x.'''];
        varStruct.(varname).string=femdatagrid;
        
    case 'FEMDATAGRIDNUM'
        varStruct.(varname).string=getbasicname('femdatagridnum');
        
    case 'FEMDATAGRIDNUMM1'
        varStruct.(varname).string=getbasicname('femdatagridnumm1');
        
    case 'FEMDATAINVMK'
        varStruct.(varname).string=getbasicname('femdatainvMK');
        
    case 'FEMDATAXCOORDM1'
        varStruct.(varname).string=getbasicname('femdataXcoordm1');
        
    case 'FEMDATAXCOORDP1'
        varStruct.(varname).string=getbasicname('femdataXcoordp1');
        
    case 'ARCVAR'
        varStruct.(varname).string=getbasicname('arcidentifiervar');
        
    case 'SUBSFLAG'
        varStruct.(varname).string='s';
        
    case 'LAGRANGEMULTCC'
        varStruct.(varname).string=getbasicname('lagrangemultcc');
        
    case 'CONTROL'
        varStruct.(varname).string=getbasicname('control');
        
    case 'SPATIALCONTROL'
        varStruct.(varname).string=['sp_' getbasicname('control')];
        
    case 'NONSPATIALCONTROL'
        varStruct.(varname).string=['nonsp_' getbasicname('control')];
        
    case 'EXISTSPATIALCONTROL'
        out=retrieveppdemodelinformation(ocStruct,'spatialcontroldependenceindex');
        varStruct.(varname).string=int2str(~isempty(out.value));
        
    case 'EXISTNONSPATIALCONTROL'
        out=retrieveppdemodelinformation(ocStruct,'nonspatialcontroldependenceindex');
        varStruct.(varname).string=int2str(~isempty(out.value));
        
    case 'LEFTINTERVALLIMIT'
        leftintervallimit=retrieveppdemodelinformation(ocStruct,'leftintervallimit');
        varStruct.(varname).string=leftintervallimit.value;
        
    case 'RIGHTINTERVALLIMIT'
        rightintervallimit=retrieveppdemodelinformation(ocStruct,'rightintervallimit');
        varStruct.(varname).string=rightintervallimit.value;

    case 'STATE'
        varStruct.(varname).string=getbasicname('state');
        
    case 'COSTATE'
        varStruct.(varname).string=getbasicname('costate');
        
    case 'EXOGENOUSFUNCTION'
        varStruct.(varname).string=getbasicname('exogenousfunction');
        
    case 'GENERATIONDATE'
        varStruct.(varname).string=datestr(now);
        
    case 'MODELNAME'
        varStruct.(varname).string=modelname(ocStruct);
        
    case 'ARCIDENTIFIER'
        varStruct.(varname).arcdependent=1;
        arcidentifier=retrieveppdemodelinformation(ocStruct,'arcidentifier');
        varStruct.(varname).string=arcidentifier.value;

    case 'ISAUTONOMOUS'
        autonomousflag=retrieveppdemodelinformation(ocStruct,'autonomous');
        varStruct.(varname).string=num2str(autonomousflag.value);

    case 'ISEXPLICITSPATIAL'
        isexplicitspatial=retrieveppdemodelinformation(ocStruct,'isexplicitspatial');
        varStruct.(varname).string=num2str(isexplicitspatial.value);
        
    case 'NUMBEROFARCS'
        arcnum=retrieveppdemodelinformation(ocStruct,'arcnum');
        varStruct.(varname).string=num2str(arcnum.value);
        
    case 'ARCARGUMENT'
        arcargument=retrieveppdemodelinformation(ocStruct,'argument');
        varStruct.(varname).string=mat2str(arcargument.value);
        
    case 'ARCCOLOR'
        arcnum=retrieveppdemodelinformation(ocStruct,'arcnum');
        colour=get(0,'DefaultAxesColorOrder');
        if arcnum.value<=size(colour,1)
            varStruct.(varname).string=mat2str(colour(1:arcnum.value,:));
        else
            varStruct.(varname).string=mat2str([colour;rand(arcnum.value-size(colour,1),3)]);
        end
        
    case 'OOPDEOBJ'
        varStruct.(varname).string='oopdeobj';
        
    case 'KRONECKERMASSMATRIX'
        kroneckermassmatrix=retrieveppdemodelinformation(ocStruct,'kroneckermassmatrix');
        varStruct.(varname).string=kroneckermassmatrix.value;

    case 'KRONECKERMASSMATRIXDEFINITION'
        kroneckermassmatrix=retrieveppdemodelinformation(ocStruct,'kroneckermassmatrix');
        pdenum=retrieveppdemodelinformation(ocStruct,'pdenum');
        varStruct.(varname).string=[kroneckermassmatrix.value '=eye(' int2str(2*pdenum.value) ');'];
        
    case 'KRONECKERDIFFUSIONMATRIX'
        kroneckermassmatrix=retrieveppdemodelinformation(ocStruct,'kroneckerdiffusionmatrix');
        varStruct.(varname).string=kroneckermassmatrix.value;
        
    case 'KRONECKERDIFFUSIONMATRIXDEFINITION'
        kroneckerdiffusionmatrix=retrieveppdemodelinformation(ocStruct,'kroneckerdiffusionmatrix');
        statedynamicsdiffusion=retrieveppdemodelinformation(ocStruct,'statedynamicsdiffusion');
        adjointsystemdiffusion=retrieveppdemodelinformation(ocStruct,'adjointsystemdiffusion');
        pdenum=retrieveppdemodelinformation(ocStruct,'pdenum');
        if length(statedynamicsdiffusion.value)==pdenum.value
            arrays='[';
            arraycs='';
            for ii=1:pdenum.value
                arrays=[arrays statedynamicsdiffusion.value{ii} ','];
                arraycs=[arraycs adjointsystemdiffusion.value{ii}.term ','];
            end
            arraycs(end)=']';
            
            varStruct.(varname).string=[kroneckerdiffusionmatrix.value '=diag(' arrays arraycs ');'];
        else
            varStruct.(varname).string='';
        end
        
    case 'EDGES'
        arcargument=retrieveppdemodelinformation(ocStruct,'argument');
        arcnum=retrieveppdemodelinformation(ocStruct,'arcnum');
        if arcnum.value>1
            nc=nchoosek(arcargument.value,2);
            varStruct.(varname).string=mat2str([nc;nc(:,[2 1])]);
        else
            varStruct.(varname).string='[]';
        end
        
    case 'CONTROLVALUESUBSTITUTION'
        controlnum=retrieveppdemodelinformation(ocStruct,'controlnum');
        controlname=retrieveppdemodelinformation(ocStruct,'controlname');
        inequalitycontrolconstraintnum=retrieveppdemodelinformation(ocStruct,'inequalitycontrolconstraintnum');
        lagrangemultipliercontrolname=retrieveppdemodelinformation(ocStruct,'lagrangemultipliercontrolname');
        indexstr='';
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
        
    case 'CONSTRAINTNUM'
        inequalitycontrolconstraintnum=retrieveppdemodelinformation(ocStruct,'inequalitycontrolconstraintnum');
        varStruct.(varname).string=num2str(inequalitycontrolconstraintnum.value);
        
    case 'CONTROLCONSTRAINTNUM'
        inequalitycontrolconstraintnum=retrieveppdemodelinformation(ocStruct,'inequalitycontrolconstraintnum');
        varStruct.(varname).string=num2str(inequalitycontrolconstraintnum.value);
        
    case 'NUMERICJACOBIAN'
        implicitcontrolsflag=implicitcontrols(ocStruct);
        varStruct.(varname).string=num2str(implicitcontrolsflag);
        
    case 'NUMERICHESSIAN'
        implicitcontrolsflag=implicitcontrols(ocStruct);
        varStruct.(varname).string=num2str(implicitcontrolsflag);
        
    case 'EXOGENOUSFUNCTIONNUM'
        exogenousfunctionnum=retrieveppdemodelinformation(ocStruct,'exogenousfunctionnum');
        varStruct.(varname).string=num2str(exogenousfunctionnum.value);
        
    case 'CONSTRAINT'
        try
            constraint=algebraicppdeterm2string(ocStruct,'constraint',1);
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
        
    case 'ISNONSPATIALCONTROL'
        nonspatialcontrol=retrieveppdemodelinformation(ocStruct,'nonspatialcontrol');
        if isempty(nonspatialcontrol.value)
            varStruct.(varname).string='0';
        else
            varStruct.(varname).string='1';
        end

    case 'STATECOSTATECOORD'
        statecostatecoordinate=retrieveppdemodelinformation(ocStruct,'statecostatecoordinate');
        varStruct.(varname).string=statecostatecoordinate.value;
        
    case 'CANONICALSYSTEMEQUATIONCOORD'
        canonicalsystemequationnum=retrieveppdemodelinformation(ocStruct,'canonicalsystemequationnum');
        varStruct.(varname).string=['1:' int2str(canonicalsystemequationnum.value)];
        
    case 'STATECOSTATENUM'
        statenum=retrieveppdemodelinformation(ocStruct,'statenum');
        costatenum=retrieveppdemodelinformation(ocStruct,'costatenum');
        varStruct.(varname).string=int2str(statenum.value+costatenum.value);
        
    case 'STATECOORD'
        statenum=retrieveppdemodelinformation(ocStruct,'statenum');
        varStruct.(varname).string=['1:' int2str(statenum.value)];
        
    case 'STATENUM'
        statenum=retrieveppdemodelinformation(ocStruct,'statenum');
        varStruct.(varname).string=int2str(statenum.value);
        
    case 'SPACEDIMENSION'
        statenum=retrieveppdemodelinformation(ocStruct,'spacedimension');
        varStruct.(varname).string=int2str(statenum.value);
        
    case 'LATEXUSERDEPENDENTNAMEFIRST'
        statename=retrieveppdemodelinformation(ocStruct,'statename');
        [val1,val2,val3]=regexp(statename.value{1},'([0-9]+)\>','split');
        if verLessThan('symbolic','8')
            if ~isempty(val2)
                varStruct.(varname).string=['$' latex(mystr2sym(val1{1})) '_' statename.value{1}(val2:val3) '$'];
            else
                varStruct.(varname).string=['$' latex(mystr2sym(val1{1})) '$'];
            end
        else
            if ~isempty(val2)
                varStruct.(varname).string=['$' latex(mystr2sym(val1{1})) '_' statename.value{1}(val2:val3) '$'];
            else
                varStruct.(varname).string=['$' latex(mystr2sym(val1{1})) '$'];
            end
        end
        %varStruct.(varname).string=['$' regexprep(statename.value{1},'([0-9]+)\>','_$1') '$'];

    case 'LATEXUSERDEPENDENTNAMESECOND'
        statenum=retrieveppdemodelinformation(ocStruct,'statenum');
        if statenum.value>1
            statename=retrieveppdemodelinformation(ocStruct,'statename');
            labelname=statename.value{2};
        else
            costatename=retrieveppdemodelinformation(ocStruct,'costatename');
            labelname=costatename.value{1};
        end
        [val1,val2,val3]=regexp(labelname,'([0-9]+)\>','split');
        if verLessThan('symbolic','8')
            if ~isempty(val2)
                varStruct.(varname).string=['$' latex(mystr2sym(val1{1})) '_' labelname(val2:val3) '$'];
            else
                varStruct.(varname).string=['$' latex(mystr2sym(val1{1})) '$'];
            end
        else
            if ~isempty(val2)
                varStruct.(varname).string=['$' latex(mystr2sym(val1{1})) '_' labelname(val2:val3) '$'];
            else
                varStruct.(varname).string=['$' latex(mystr2sym(val1{1})) '$'];
            end
        end
        
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
        
    case 'COSTATECOORD'
        statenum=retrieveppdemodelinformation(ocStruct,'statenum');
        costatenum=retrieveppdemodelinformation(ocStruct,'costatenum');
        varStruct.(varname).string=[int2str(statenum.value+1) ':' int2str(statenum.value+costatenum.value)];
        
    case 'ODEDIM'
        varStruct.(varname).arcdependent=1;
        odedim=retrieveppdemodelinformation(ocStruct,'odedim');
        arcnum=retrieveppdemodelinformation(ocStruct,'arcnum');
        for ii=1:arcnum.value
            varStruct.(varname).string{ii}=int2str(odedim.value);
        end
        
        
    case 'PARAMETERVALUES'
        % the strings used for assigning parametervalues from the par
        % argument in functions to the spcific paraemter variable names,
        % e.g. r=pararg(1); ...
        varStruct.(varname).multline=2;
        basename=getbasicname('parametervariables');
        parvar=algebraicppdeterm2string(ocStruct,'parametervariables');
        for ii=1:numel(parvar.term)
            varStruct.(varname).string{ii}=[parvar.term{ii} '=' basename '(' num2str(ii) ');'];
        end
        
    case 'EXOGENOUSFUNCTIONTERM'
        try
            exogenousfunctionterm=algebraicppdeterm2string(ocStruct,'exogenousfunctionterm',1,'');
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
            canonicalsystem=algebraicppdeterm2string(ocStruct,'canonicalsystem',1,'');
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
        
    case 'SYMBOLICCANONICALSYSTEMDYNAMICS'
        try
            canonicalsystem=algebraicppdeterm2string(ocStruct,'symboliccanonicalsystem',0,'');
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
            arcnum=retrieveppdemodelinformation(ocStruct,'arcnum');
            varStruct.(varname).arcidentifier=field2arcidentifier(arcfield);
            for ii=1:numel(arcfield)
                equilibriumequation=algebraicppdeterm2string(ocStruct,'equilibriumequation',0,arcfield{ii});
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
        
    case 'CANONICALSYSTEMJACOBIAN'
        try
            arcfield=getarcclass(ocStruct,'canonicalsystemjacobian');
            arcnum=retrieveppdemodelinformation(ocStruct,'arcnum');
            varStruct.(varname).arcidentifier=field2arcidentifier(arcfield);
            for ii=1:numel(arcfield)
                canonicalsystemjacobian=algebraicppdeterm2string(ocStruct,'canonicalsystemjacobian',2,arcfield{ii});
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
    case 'SPECIFICPONTRYAGINFUNCTION'
        try
            arcfield=getarcclass(ocStruct,'canonicalsystemjacobian');
            arcnum=retrieveppdemodelinformation(ocStruct,'arcnum');
            varStruct.(varname).arcidentifier=field2arcidentifier(arcfield);
            for ii=1:numel(arcfield)
                specificpontryaginfunction=algebraicppdeterm2string(ocStruct,'specificpontryaginfunction',1,arcfield{ii});
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
            arcnum=retrieveppdemodelinformation(ocStruct,'arcnum');
            varStruct.(varname).arcidentifier=field2arcidentifier(arcfield);
            for ii=1:numel(arcfield)
                canonicalsystemjacobian=algebraicppdeterm2string(ocStruct,'symboliccanonicalsystemjacobian',0,arcfield{ii});
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
            arcnum=retrieveppdemodelinformation(ocStruct,'arcnum');
            varStruct.(varname).arcidentifier=field2arcidentifier(arcfield);
            for ii=1:numel(arcfield)
                canonicalsystemjacobian=algebraicppdeterm2string(ocStruct,'symbolicstatecontrolcanonicalsystemjacobian',0,arcfield{ii});
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
            arcnum=retrieveppdemodelinformation(ocStruct,'arcnum');
            varStruct.(varname).arcidentifier=field2arcidentifier(arcfield);
            for ii=1:numel(arcfield)
                canonicalsystemjacobian=algebraicppdeterm2string(ocStruct,'symbolicstatecontrolcanonicalsystem',0,arcfield{ii});
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
                canonicalsystemparameterjacobian=algebraicppdeterm2string(ocStruct,'canonicalsystemparameterjacobian',2,arcfield{ii});
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
            arcnum=retrieveppdemodelinformation(ocStruct,'arcnum');
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
                canonicalsystemparameterjacobian=algebraicppdeterm2string(ocStruct,'statecontrolcanonicalsystemparameterjacobian',0,arcfield{ii});
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
            arcnum=retrieveppdemodelinformation(ocStruct,'arcnum');
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
                canonicalsystemderivativetime=algebraicppdeterm2string(ocStruct,'canonicalsystemderivativetime',0,arcfield{ii});
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
            arcnum=retrieveppdemodelinformation(ocStruct,'arcnum');
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
                canonicalsystemhessian=algebraicppdeterm2string(ocStruct,'canonicalsystemhessian',0,arcfield{ii});
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
            arcnum=retrieveppdemodelinformation(ocStruct,'arcnum');
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
                canonicalsystemparameterhessian=algebraicppdeterm2string(ocStruct,'canonicalsystemtotalhessian',0,arcfield{ii});
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
            arcnum=retrieveppdemodelinformation(ocStruct,'arcnum');
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
                pontryaginfunctionDu2=algebraicppdeterm2string(ocStruct,'pontryaginfunctionDu2',0,arcfield{ii});
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
            arcnum=retrieveppdemodelinformation(ocStruct,'arcnum');
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
        parameternum=retrieveppdemodelinformation(ocStruct,'parameternum');
        for ii=1:numel(arcfield)
            equationnum=retrieveppdemodelinformation(ocStruct,'canonicalsystemequationnum',field2arcidentifier(arcfield{ii}));
            varStruct.(varname).string{ii}=['1:' int2str(equationnum.value) ',1:'  int2str(equationnum.value) ',' int2str(equationnum.value+1) ':' int2str(parameternum.value+equationnum.value)];
            arcnum=retrieveppdemodelinformation(ocStruct,'arcnum');
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
                spatialcontrolvalue=algebraicppdeterm2string(ocStruct,'spatialcontrolvalue',1,arcfield{ii});
                for jj=1:numel(spatialcontrolvalue.term)
                    if jj==1
                        varStruct.(varname).string{ii}{jj}=['out=' spatialcontrolvalue.term{jj}];
                    else
                        varStruct.(varname).string{ii}{jj}=spatialcontrolvalue.term{jj};
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
        
    case 'NONSPATIALOPTIMALCONTROLVALUE'
        try
            arcfield=getarcclass(ocStruct,'controlvalue');
            varStruct.(varname).arcidentifier=field2arcidentifier(arcfield);
            for ii=1:numel(arcfield)
                nonspatialcontrolvalue=algebraicppdeterm2string(ocStruct,'nonspatialcontrolvalue',1,arcfield{ii});
                for jj=1:numel(nonspatialcontrolvalue.term)
                    if jj==1
                        varStruct.(varname).string{ii}{jj}=['out=' nonspatialcontrolvalue.term{jj}];
                    else
                        varStruct.(varname).string{ii}{jj}=nonspatialcontrolvalue.term{jj};
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
        
    case 'OPTIMALCONTROLVALUEIMPLICIT'
        try
            arcfield=getarcclass(ocStruct,'controlvalue');
            varStruct.(varname).arcidentifier=field2arcidentifier(arcfield);
            for ii=1:numel(arcfield)
                controlvalue=algebraicppdeterm2string(ocStruct,'controlvalue',1,arcfield{ii});
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
                controlvalue=algebraicppdeterm2string(ocStruct,'optimalvalue4statecontrol',1,arcfield{ii});
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
            arcnum=retrieveppdemodelinformation(ocStruct,'arcnum');
            varStruct.(varname).arcidentifier=field2arcidentifier(arcfield);
            for ii=1:numel(arcfield)
                optimalcostatevalue=algebraicppdeterm2string(ocStruct,'optimalcostatevalue',1,arcfield{ii});
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
            arcnum=retrieveppdemodelinformation(ocStruct,'arcnum');
            varStruct.(varname).arcidentifier=field2arcidentifier(arcfield);
            for ii=1:numel(arcfield)
                explicitstatevalue=algebraicppdeterm2string(ocStruct,'explicitstatevalue',1,arcfield{ii});
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
                symbolicstatevalue=algebraicppdeterm2string(ocStruct,'symbolicstatevalue',0,arcfield{ii});
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
                controlvalue=algebraicppdeterm2string(ocStruct,'symboliccontrolvalue',0,arcfield{ii});
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
        
    case 'SYMBOLICOPTIMALCONTROLVALUE4STATECONTROL'
        try
            arcfield=getarcclass(ocStruct,'controlvalue');
            varStruct.(varname).arcidentifier=field2arcidentifier(arcfield);
            for ii=1:numel(arcfield)
                controlvalue=algebraicppdeterm2string(ocStruct,'symboliccontrolvalue4statecontrol',0,arcfield{ii});
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
                controlvalue=algebraicppdeterm2string(ocStruct,'symboliccostatevalue',0,arcfield{ii});
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
            if isempty(arcid)
                varStruct.(varname).string{1}='';
            end
            for ii=1:numel(arcid)
                lagrangemultcc=algebraicppdeterm2string(ocStruct,'lagrangemultcc',1,arcid{ii});
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
                    lagrangemultsc=algebraicppdeterm2string(ocStruct,'lagrangemultsc',1,arcid{ii});
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
                lagrangemultcc=algebraicppdeterm2string(ocStruct,'symboliclagrangemultcc',0,arcid{ii});
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
                lagrangemultsc=algebraicppdeterm2string(ocStruct,'symboliclagrangemultsc',0,arcid{ii});
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
        
    case 'SYMBOLICDISCOBJECTIVEFUNCTION'
        try
            objectivefunction=algebraicppdeterm2string(ocStruct,'symbolicdiscobjectivefunction',0);
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
            objectivefunction=algebraicppdeterm2string(ocStruct,'discobjectivefunction',1);
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
                discobjectivefunctionderivativetime=algebraicppdeterm2string(ocStruct,'discobjectivefunctionderivativetime',0,arcfield{ii});
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
            arcnum=retrieveppdemodelinformation(ocStruct,'arcnum');
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
            arcnum=retrieveppdemodelinformation(ocStruct,'arcnum');
            varStruct.(varname).arcidentifier=field2arcidentifier(arcfield);
            for ii=1:arcnum.value
                discobjectivefunctionDx=algebraicppdeterm2string(ocStruct,'discobjectivefunctionjacobian',2,arcfield{ii});
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
            arcnum=retrieveppdemodelinformation(ocStruct,'arcnum');
            varStruct.(varname).arcidentifier=field2arcidentifier(arcfield);
            for ii=1:arcnum.value
                discobjectivefunctionDx=algebraicppdeterm2string(ocStruct,'discobjectivefunctionparameterjacobian',2,arcfield{ii});
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
        
    case 'HAMILTONIANFUNCTION'
        try
            hamiltonianfunction=algebraicppdeterm2string(ocStruct,'hamiltonianfunction',1);
            varStruct.(varname).string{1}=['out=' hamiltonianfunction.term{1} ';'];
            varStruct.(varname).multline=1;
            varStruct.(varname).vectorize=1;
            varStruct.(varname).type='algterm';
        catch
            %lasterr
            varStruct.(varname).string=varname;
        end
        
    case 'PONTRYAGINFUNCTION'
        try
            pontryaginfunction=algebraicppdeterm2string(ocStruct,'pontryaginfunction',1);
            varStruct.(varname).string{1}=['out=' pontryaginfunction.term{1} ';'];
            varStruct.(varname).multline=1;
            varStruct.(varname).vectorize=1;
            varStruct.(varname).type='algterm';
        catch
            %lasterr
            varStruct.(varname).string=varname;
        end
        
    case 'SALVAGEVALUE'
        try
            salvagevalue=algebraicppdeterm2string(ocStruct,'salvagevalue',1);
            varStruct.(varname).string{1}=['out=' salvagevalue.term ';'];
            varStruct.(varname).multline=0;
            varStruct.(varname).vectorize=0;
            varStruct.(varname).type='algterm';
        catch
            %lasterr
            varStruct.(varname).string=varname;
        end
        
    case 'DISCOUNTRATE'
        try
            discountratevariable=retrieveppdemodelinformation(ocStruct,'discountratevariable');
            varStruct.(varname).string{1}=discountratevariable.value;
            varStruct.(varname).multline=0;
            varStruct.(varname).vectorize=0;
            varStruct.(varname).type='algterm';
        catch
            %lasterr
            varStruct.(varname).string=varname;
        end
        
    case 'DISCOUNTEDSALVAGEVALUE'
        try
            salvagevalue=algebraicppdeterm2string(ocStruct,'discountedsalvagevalue',1);
            varStruct.(varname).string{1}=['out=' salvagevalue.term ';'];
            varStruct.(varname).multline=0;
            varStruct.(varname).vectorize=0;
            varStruct.(varname).type='algterm';
        catch
            %lasterr
            varStruct.(varname).string=varname;
        end
        
    case 'SYMBOLICOBJECTIVEFUNCTION'
        try
            symbolicobjectivefunction=algebraicppdeterm2string(ocStruct,'symbolicobjectivefunction',0);
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
            symbolicpontryaginfunction=algebraicppdeterm2string(ocStruct,'symbolicpontryaginfunction',0);
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
            arcnum=retrieveppdemodelinformation(ocStruct,'arcnum');
            varStruct.(varname).arcidentifier=field2arcidentifier(arcfield);
            for ii=1:arcnum.value
                pontryaginfunctionDx=algebraicppdeterm2string(ocStruct,'pontryaginfunctionDX',0,arcfield{ii});
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
            hamiltonianDu=algebraicppdeterm2string(ocStruct,'hamiltonianDu',1);
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
            statedynamics=algebraicppdeterm2string(ocStruct,'statedynamics',1,'');
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
            adjointsystem=algebraicppdeterm2string(ocStruct,'adjointsystem',1,'');
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
        
    case 'OBJECTIVEFUNCTION'
        try
            objectivefunction=algebraicppdeterm2string(ocStruct,'objectivefunction',1);
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
            transversalitycondition=algebraicppdeterm2string(ocStruct,'transversalitycondition',0);
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