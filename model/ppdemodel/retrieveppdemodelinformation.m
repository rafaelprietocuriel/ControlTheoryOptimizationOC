function data=retrieveppdemodelinformation(ocStruct,propertyclass,arcidentifier,symkernel)
%
% RETRIEVEMODELINFORMATION returns a structure with information about the
% ocmat model
% this function is the interface to the structure OCSTRUCT as derived
% from the models initialization file. Its purpose is to allow the same
% commands even if the structure OCSTRUCT is changed. In that case only
% 'retrieveodemodelinformation' has to be adapted.

data.type='';
data.value='';
data.description='';

if isempty(ocStruct)
    return
end
if nargin==2
    arcidentifier='';
end

switch propertyclass
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % GENERAL
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    case 'arcidentifier'
        % arc identifiers specify different combinations of (in)active
        % constraints and/or multiple solutions of the  Hamiltonian
        % maximizing condition
        data.type='char';
        data.value=ocStruct.arc.identifier;
        data.description='identifier to differentiate between specifications of the canonical system';

    case 'argument'
        % arcidentifiers (char) are transformed the corresponding arguments
        % (numeric)
        data.type='integervector';
        data.value=ocStruct.arc.argument;
        data.description='arcarguments is a vector containing the numeric values of the arcidentifiers';

    case 'independent'
        % variable name of the indepent variable (usually time)
        data.type='char';
        data.value=ocStruct.variable.independent.name;
        data.description='independent variable';

    case 'time'
        % variable name of the indepent variable (usually time)
        data.type='char';
        data.value=ocStruct.variable.independent.time.name;
        data.description='independent time variable';

    case 'space'
        % variable name of the indepent variable (usually time)
        data.type='char';
        data.value=ocStruct.variable.independent.space.name;
        data.description='independent space variable';

    case 'spacemid'
        % variable name of the indepent variable (usually time)
        data.type='char';
        data.value=[ocStruct.variable.independent.space.name 'mid'];
        data.description='independent space variable';

    case 'autonomous'
        % variable name of the indepent variable (usually time)
        data.type='boolean';
        data.value=ocStruct.variable.independent.property.autonomous;
        data.description='model not explicitly depending on time';

    case 'isexplicitspatial'
        % variable name of the indepent variable (usually time)
        data.type='boolean';
        data.value=ocStruct.variable.independent.property.explicitspatial;
        data.description='model (not) depending explicitly on space';
        
    case 'femdataucoord'
        % coordinates of the spatial and nonspatial state and costate
        controlnum=retrieveppdemodelinformation(ocStruct,'controlnum');
        gridnum=getbasicname('femdatagridnum');
        spatialcontrolpendenceindex=retrieveppdemodelinformation(ocStruct,'spatialcontroldependenceindex');
        nonspatialcontroldependenceindex=retrieveppdemodelinformation(ocStruct,'nonspatialcontroldependenceindex');
        
        value=cell(1,controlnum.value);
        ctr=0;
        base='';
        for ii=1:controlnum.value
            if any(spatialcontrolpendenceindex.value==ii)
                if ii==1
                    value{ii}=['1:' gridnum];
                elseif ii==2
                    value{ii}=[gridnum '+1:2*' gridnum];
                else
                    value{ii}=[int2str(ii-1) '*' gridnum '+1:' int2str(ii) '*' gridnum];
                end
            elseif any(nonspatialcontroldependenceindex.value==ii)
                if ctr==0 && ii>1
                    if ii==2
                        base=gridnum;
                    else
                        base=[int2str(ii-1) '*' gridnum];
                    end
                end
                ctr=ctr+1;
                if ii==1
                    value{ii}=int2str(ctr);
                else
                value{ii}=[base '+' int2str(ctr)];
                end
            end
        end
        
        data.type='char cell';
        data.value=value;
        data.description='fem discretization index of control variables';

    case 'statecostatecoordinate'
        % coordinates of the spatial and nonspatial state and costate
        gridnum=getbasicname('femdatagridnum');
        spatialstatedependenceindex=retrieveppdemodelinformation(ocStruct,'spatialstatedependenceindex');
        nonspatialstatedependenceindex=retrieveppdemodelinformation(ocStruct,'nonspatialstatedependenceindex');
        spatialcostatedependenceindex=retrieveppdemodelinformation(ocStruct,'spatialcostatedependenceindex');
        nonspatialcostatedependenceindex=retrieveppdemodelinformation(ocStruct,'nonspatialcostatedependenceindex');
        spatialcoordnum=spatialstatedependenceindex.value(end)+spatialcostatedependenceindex.value(end);
        nonspatialcoordnum=nonspatialstatedependenceindex.value+nonspatialcostatedependenceindex.value;
        
        if isempty(nonspatialcoordnum) || isempty(nonspatialcoordnum)
            coord=['1:' int2str(spatialcoordnum) '*' gridnum];
        else
            coord=['1:' int2str(spatialcoordnum) '*' gridnum '+' int2str(nonspatialcoordnum)];
        end
        
        data.type='char';
        data.value=coord;
        data.description='fem discretization index of dependent variables';
        
    case 'femdataXcoord'
        % coordinates of the spatial and nonspatial state and costate
        statenum=retrieveppdemodelinformation(ocStruct,'statenum');
        costatenum=retrieveppdemodelinformation(ocStruct,'costatenum');
        gridnum=getbasicname('femdatagridnum');
        spatialstatedependenceindex=retrieveppdemodelinformation(ocStruct,'spatialstatedependenceindex');
        nonspatialstatedependenceindex=retrieveppdemodelinformation(ocStruct,'nonspatialstatedependenceindex');
        spatialcostatedependenceindex=retrieveppdemodelinformation(ocStruct,'spatialcostatedependenceindex');
        nonspatialcostatedependenceindex=retrieveppdemodelinformation(ocStruct,'nonspatialcostatedependenceindex');
        
        value=cell(1,statenum.value+costatenum.value);
        ctr=0;
        base='';
        for ii=1:statenum.value
            if any(spatialstatedependenceindex.value==ii)
                if ii==1
                    value{ii}=['1:' gridnum];
                elseif ii==2
                    value{ii}=[gridnum '+1:2*' gridnum];
                else
                    value{ii}=[int2str(ii-1) '*' gridnum '+1:' int2str(ii) '*' gridnum];
                end
            elseif any(nonspatialstatedependenceindex.value==ii)
                if ctr==0
                    base=[int2str(ii-1) '*' gridnum];
                end
                ctr=ctr+1;
                
                value{ii}=[base '+' int2str(ctr)];
            end
        end
        if ctr==0
            if ii==1
                base=gridnum;
            else
                base=[int2str(ii) '*' gridnum];
            end
        else
            base=[base '+' int2str(ctr)];
        end
        ctr=0;
        for ii=1:costatenum.value
            if any(spatialcostatedependenceindex.value==ii)
                if ii==1
                    value{statenum.value+ii}=[base '+(1:' gridnum ')'];
                elseif ii==2
                    value{statenum.value+ii}=[base '+(' gridnum '+1:2*' gridnum ')'];
                else
                    value{statenum.value+ii}=[base '+(' int2str(ii-1) '*' gridnum '+1:' int2str(ii) '*' gridnum ')'];
                end
            elseif any(nonspatialcostatedependenceindex.value==ii)
                ctr=ctr+1;
                if ctr==1
                    base=[base '+' int2str(ii-1) '*' gridnum];
                end
                value{statenum.value+ii}=[base '+' int2str(ctr)];
            end
        end
        
        data.type='char cell';
        data.value=value;
        data.description='fem discretization index of dependent variables';
        
    case 'femdataXcoordm1'
        % coordinates of the spatial and nonspatial state and costate
        statenum=retrieveppdemodelinformation(ocStruct,'statenum');
        costatenum=retrieveppdemodelinformation(ocStruct,'costatenum');
        gridnum=getbasicname('femdatagridnum');
        spatialstatedependenceindex=retrieveppdemodelinformation(ocStruct,'spatialstatedependenceindex');
        nonspatialstatedependenceindex=retrieveppdemodelinformation(ocStruct,'nonspatialstatedependenceindex');
        spatialcostatedependenceindex=retrieveppdemodelinformation(ocStruct,'spatialcostatedependenceindex');
        nonspatialcostatedependenceindex=retrieveppdemodelinformation(ocStruct,'nonspatialcostatedependenceindex');
        
        value=cell(1,statenum.value+costatenum.value);
        ctr=0;
        base='';
        for ii=1:statenum.value
            if any(spatialstatedependenceindex.value==ii)
                if ii==1
                    value{ii}=['1:' gridnum '-1'];
                elseif ii==2
                    value{ii}=[gridnum '+1:2*' gridnum '-1'];
                else
                    value{ii}=[int2str(ii-1) '*' gridnum '+1:' int2str(ii) '*' gridnum '-1'];
                end
            elseif any(nonspatialstatedependenceindex.value==ii)
                if ctr==0
                    base=[int2str(ii-1) '*' gridnum];
                end
                ctr=ctr+1;
                
                value{ii}=[base '+(' int2str(ctr) ':' base '+' int2str(ctr) ')'];
            end
        end
        if ctr==0
            if ii==1
                base=gridnum;
            else
                base=[int2str(ii) '*' gridnum];
            end
        else
            base=[base '+' int2str(ctr)];
        end
        ctr=0;
        for ii=1:costatenum.value
            if any(spatialcostatedependenceindex.value==ii)
                if ii==1
                    value{statenum.value+ii}=[base '+(1:' gridnum '-1)'];
                elseif ii==2
                    value{statenum.value+ii}=[base '+(' gridnum '+1:2*' gridnum '-1)'];
                else
                    value{statenum.value+ii}=[base '+(' int2str(ii-1) '*' gridnum '+1:' int2str(ii) '*' gridnum '-1)'];
                end
            elseif any(nonspatialcostatedependenceindex.value==ii)
                ctr=ctr+1;
                if ctr==1
                    base=[base '+' int2str(ii-1) '*' gridnum];
                end
                value{statenum.value+ii}=[base '+(' int2str(ctr) ':' base '+' int2str(ctr) ')'];
            end
        end
        
        data.type='char cell';
        data.value=value;
        data.description='fem discretization index of dependent variables';
        
    case 'femdataXcoordp1'
        % coordinates of the spatial and nonspatial state and costate
        statenum=retrieveppdemodelinformation(ocStruct,'statenum');
        costatenum=retrieveppdemodelinformation(ocStruct,'costatenum');
        gridnum=getbasicname('femdatagridnum');
        spatialstatedependenceindex=retrieveppdemodelinformation(ocStruct,'spatialstatedependenceindex');
        nonspatialstatedependenceindex=retrieveppdemodelinformation(ocStruct,'nonspatialstatedependenceindex');
        spatialcostatedependenceindex=retrieveppdemodelinformation(ocStruct,'spatialcostatedependenceindex');
        nonspatialcostatedependenceindex=retrieveppdemodelinformation(ocStruct,'nonspatialcostatedependenceindex');
        
        value=cell(1,statenum.value+costatenum.value);
        ctr=0;
        base='';
        for ii=1:statenum.value
            if any(spatialstatedependenceindex.value==ii)
                if ii==1
                    value{ii}=['2:' gridnum];
                elseif ii==2
                    value{ii}=[gridnum '+2:2*' gridnum];
                else
                    value{ii}=[int2str(ii-1) '*' gridnum '+2:' int2str(ii) '*' gridnum];
                end
            elseif any(nonspatialstatedependenceindex.value==ii)
                if ctr==0
                    base=[int2str(ii-1) '*' gridnum];
                end
                ctr=ctr+1;
                
                value{ii}=[base '+(' int2str(ctr) ':' base '+' int2str(ctr) ')'];
            end
        end
        if ctr==0
            if ii==1
                base=gridnum;
            else
                base=[int2str(ii) '*' gridnum];
            end
        else
            base=[base '+' int2str(ctr)];
        end
        ctr=0;
        for ii=1:costatenum.value
            if any(spatialcostatedependenceindex.value==ii)
                if ii==1
                    value{statenum.value+ii}=[base '+(2:' gridnum ')'];
                elseif ii==2
                    value{statenum.value+ii}=[base '+(' gridnum '+2:2*' gridnum ')'];
                else
                    value{statenum.value+ii}=[base '+(' int2str(ii-1) '*' gridnum '+2:' int2str(ii) '*' gridnum ')'];
                end
            elseif any(nonspatialcostatedependenceindex.value==ii)
                ctr=ctr+1;
                if ctr==1
                    base=[base '+' int2str(ii-1) '*' gridnum];
                end
                value{statenum.value+ii}=[base '+(' int2str(ctr) ':' base '+' int2str(ctr) ')'];
            end
        end
        
        data.type='char cell';
        data.value=value;
        data.description='fem discretization index of dependent variables';
        
    case 'kroneckermassmatrix'
        % variable name of the state(s) as defined by the user in the
        % initialization file
        data.type='char';
        data.value='kronMatMass';
        data.description='variable of the matrix for the Kroneckerproduct of the mass matrix';
        
    case 'spacedimension'
        % variable name of the state(s) as defined by the user in the
        % initialization file
        data.type='integer';
        if isfield(ocStruct.spacegeometry,'R1') && isfield(ocStruct.spacegeometry.R1,'interval')
            data.value=1;
        else
            data.value='';
        end
        data.description='dimension of the space';
        
    case 'kroneckerdiffusionmatrix'
        % variable name of the state(s) as defined by the user in the
        % initialization file
        data.type='char';
        data.value='kronMatDiffusion';
        data.description='variable of the matrix for the Kroneckerproduct of the mass matrix';

    case 'leftintervallimit'
        % dependent variable names of the Hamiltonian maximizing condition
        if isfield(ocStruct.spacegeometry,'R1') && isfield(ocStruct.spacegeometry.R1,'interval')
            data.value=ocStruct.spacegeometry.R1.interval.boundary.left.term;
        else
            data.value='';
        end
        data.type='char';
        data.description='if the geometry of the space is an interval the left limit is returned.';

    case 'rightintervallimit'
        % dependent variable names of the Hamiltonian maximizing condition
        if isfield(ocStruct.spacegeometry,'R1') && isfield(ocStruct.spacegeometry.R1,'interval')
        data.value=ocStruct.spacegeometry.R1.interval.boundary.right.term;
        else
            data.value='';
        end
        data.type='char';
        data.description='if the geometry of the space is an interval the left limit is returned.';
        
    case 'statename'
        % variable name of the state(s) as defined by the user in the
        % initialization file
        data.type='cellchar';
        data.value={ocStruct.variable.state.name};
        data.description='name of state';

    case 'statedependence'
        % variable name of the state(s) as defined by the user in the
        % initialization file
        data.type='cellchar';
        data.value={ocStruct.variable.state.dependence};
        data.description='time/space dependence of state';

    case 'spatialstatedependenceindex'
        % variable name of the state(s) as defined by the user in the
        % initialization file
        statedependence=retrieveppdemodelinformation(ocStruct,'statedependence');
        data.type='integer array';
        data.value=find(strcmp(statedependence.value,'2'));
        data.description='index of spatially dependent state';

    case 'nonspatialstatedependenceindex'
        % variable name of the state(s) as defined by the user in the
        % initialization file
        statedependence=retrieveppdemodelinformation(ocStruct,'statedependence');
        data.type='integer array';
        data.value=find(strcmp(statedependence.value,'0'));
        data.description='index of spatially independent state';

    case 'spatialcontroldependenceindex'
        % variable name of the state(s) as defined by the user in the
        % initialization file
        controldependence=retrieveppdemodelinformation(ocStruct,'controldependence');
        data.type='integer array';
        data.value=find(strcmp(controldependence.value,'2'));
        data.description='index of spatially dependent control';

    case 'nonspatialcontroldependenceindex'
        % variable name of the state(s) as defined by the user in the
        % initialization file
        controldependence=retrieveppdemodelinformation(ocStruct,'controldependence');
        data.type='integer array';
        data.value=find(strcmp(controldependence.value,'0'));
        data.description='index of spatially independent control';

    case 'spatialcostatedependenceindex'
        % variable name of the state(s) as defined by the user in the
        % initialization file
        costatedependence=retrieveppdemodelinformation(ocStruct,'costatedependence');
        data.type='integer array';
        data.value=find(strcmp(costatedependence.value,'2'));
        data.description='index of spatially dependent state';

    case 'nonspatialcostatedependenceindex'
        % variable name of the state(s) as defined by the user in the
        % initialization file
        costatedependence=retrieveppdemodelinformation(ocStruct,'costatedependence');
        data.type='integer array';
        data.value=find(strcmp(costatedependence.value,'0'));
        data.description='index of spatially independent state';

    case 'stateboundarycondition'
        % variable name of the state(s) as defined by the user in the
        % initialization file
        data.type='cellchar';
        data.value={ocStruct.constraint.derivative.state.boundarycondition};
        data.description='spatial boundary conditions of the state(s)';

    case 'controldependence'
        % variable name of the state(s) as defined by the user in the
        % initialization file
        data.type='cellchar';
        data.value={ocStruct.variable.control.dependence};
        data.description='time/space dependence of control';

    case 'nonspatialcontrol'
        % index of controls that do not explicitly depende on spatial
        % dimension
        controldependence=retrieveppdemodelinformation(ocStruct,'controldependence');
        data.type='integer array';
        data.value=find(strcmp(controldependence.value,'0'));
        data.description='index of nonspatial control';
        
    case 'costatedependence'
        % variable name of the state(s) as defined by the user in the
        % initialization file
        data.type='cellchar';
        data.value={ocStruct.variable.costate.dependence};
        data.description='time/space dependence of costate';

    case 'statedynamics'
        % ODEs describing the evolution of the state(s)
        data.type='mathchar';
        data.value={ocStruct.constraint.derivative.state.term};
        data.description='ODEs/PDEs of the state dynamics';

    case 'statedynamicstype'
        % ODEs describing the evolution of the state(s)
        data.type='mathchar';
        data.value={ocStruct.constraint.derivative.state.type};
        data.description='Type of the state dynamics';

    case 'pdenum'
        % ODEs describing the evolution of the state(s)
        data.type='integer';
        statedynamicstype=retrieveppdemodelinformation(ocStruct,'statedynamicstype');
        data.value=sum(strcmp(statedynamicstype.value,'pde'));
        data.description='Number of PDEs in the state dynamics';

    case 'odenum'
        % ODEs describing the evolution of the state(s)
        data.type='integer';
        statedynamicstype=retrieveppdemodelinformation(ocStruct,'statedynamicstype');
        data.value=sum(strcmp(statedynamicstype.value,'ode'));
        data.description='Number of ODEs in the state dynamics';

    case 'statedynamicsconvection'
        % ODEs describing the evolution of the state(s)
        data.type='mathchar';
        try
            for ii=1:length(ocStruct.constraint.derivative.state)
                data.value{ii}=ocStruct.constraint.derivative.state(ii).convection.term;
            end
        catch
            data.value=[];
        end
        data.description='Convection term of the state dynamics';

    case 'statedynamicsconvectionvariable'
        % ODEs describing the evolution of the state(s)
        data.type='mathchar';
        try
            for ii=1:length(ocStruct.constraint.derivative.state)
                data.value{ii}=ocStruct.constraint.derivative.state(ii).convection.variable;
            end
        catch
            data.value=[];
        end
        data.description='Diffusion term of the state dynamics';

    case 'statedynamicsdiffusion'
        % ODEs describing the evolution of the state(s)
        data.type='mathchar';
        try
            for ii=1:length(ocStruct.constraint.derivative.state)
                data.value{ii}=ocStruct.constraint.derivative.state(ii).diffusion.term;
            end
        catch
            data.value=[];
        end
        data.description='Diffusion term of the state dynamics';

    case 'statedynamicsdiffusionvariable'
        % ODEs describing the evolution of the state(s)
        data.type='mathchar';
        try
            for ii=1:length(ocStruct.constraint.derivative.state)
                data.value{ii}=ocStruct.constraint.derivative.state(ii).diffusion.variable;
            end
        catch
            data.value=[];
        end
        data.description='Diffusion term of the state dynamics';

    case 'costatename'
        % variable name of the state(s) as defined by the user in the
        % initialization file
        data.type='cellchar';
        data.value={ocStruct.variable.costate.name};
        data.description='name of costate';
        
    case 'controlname'
        % variable name of the control(s) as defined by the user in the
        % initialization file
        data.type='cellchar';
        data.value={ocStruct.variable.control.name};
        data.description='name of control';

    case 'costatenum'
        % the number of costate(s)
        data.type='integer';
        data.value=length(ocStruct.variable.costate);
        data.description='number of costate';

    case 'parameternum'
        % number of parameter(s)
        data.type='integer';
        if isfield(ocStruct,'parameter')
            data.value=ocStruct.parameter.num;
        else
            data.value=0;
        end
        data.description='number of exogenous parameters';

    case 'objectiveintegrand'
        % integrand for objective value function without the possible
        % discounting factor
        data.type='mathchar';
        data.value=ocStruct.objective.integral.function.term;
        data.description='integrand for objective value function';

    case 'objectiveintegranddependence'
        % variable name of the state(s) as defined by the user in the
        % initialization file
        data.type='char';
        data.value=ocStruct.objective.integral.function.dependence;
        data.description='time/space dependence of objective integrand';

    case 'exogenousfunctionnum'
        % variable name of the control(s) as defined by the user in the
        % initialization file
        data.type='cellchar';
        if isfield(ocStruct,'exogenousfunction')
            data.value=length(fieldnames(ocStruct.exogenousfunction));
        else
            data.value=0;
        end
        data.description='names of exogenous functions';

    case 'exogenousfunctionname'
        % variable name of the control(s) as defined by the user in the
        % initialization file
        data.type='cellchar';
        if isfield(ocStruct,'exogenousfunction')
            value=fieldnames(ocStruct.exogenousfunction).';
            for ii=1:length(value)
                value{ii}=[value{ii} '()'];
            end
            data.value=value;
        else
            data.value=[];
        end
        data.description='names of exogenous functions';

    case 'exogenousfunctionterm'
        % variable name of the control(s) as defined by the user in the
        % initialization file
        data.type='cellchar';
        exogenousfunctionname=retrieveodemodelinformation(ocStruct,'exogenousfunctionname');
        if ~isempty(exogenousfunctionname.value)
            for ii=1:length(exogenousfunctionname.value)
                data.value{ii}=ocStruct.exogenousfunction.(exogenousfunctionname.value{ii}).term;
            end
        else
            data.value=[];
        end
        data.description='terms of exogenous functions';
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%
        % Objective Function
        %%%%%%%%%%%%%%%%%%%%%%%%%%
        

    case 'objectivetype'
        % type of the terms which build up the objective value function
        % (usually integral and Salvage value)
        data.type='mathchar';
        data.value=fieldnames(ocStruct.objective);
        data.description='type of terms which build up the objective value function';
        
    case 'discountratevariable'
        % discount variable used in the discount factor of the objective
        % value function
        if isfield(ocStruct.objective,'integral')
            data.type='char';
            data.value=ocStruct.objective.integral.discountrate;
            data.description='discount variable of the objective value function';
        end

    case 'discobjectivefunction'
        objectivefunction=retrieveppdemodelinformation(ocStruct,'objectiveintegrand');
        discountfactor=retrieveppdemodelinformation(ocStruct,'discountfactor');
        if isfield(ocStruct.objective,'integral')
            data.type='mathchar';
            data.value=[discountfactor.value '*(' objectivefunction.value ')'];
            data.description='the discounted objective integrand';
        end

    case 'specificdiscobjectivefunction'
        discobjectivefunction=retrieveppdemodelinformation(ocStruct,'discobjectivefunction');
        optimalvalue=retrieveppdemodelinformation(ocStruct,'optimalvalue',arcidentifier);
        maximizingvariable=retrieveppdemodelinformation(ocStruct,'maximizingvariable');
        for ii=1:numel(maximizingvariable.value)
            discobjectivefunction.value=ocmatsubs(discobjectivefunction.value,[maximizingvariable.value{ii} '=' optimalvalue.value.(maximizingvariable.value{ii})],symkernel);
        end
        data.type='mathchar';
        data.value=discobjectivefunction.value;
        data.description='the discounted objective integrand';

    case 'discobjectivefunctionDX'
        costatename=retrieveppdemodelinformation(ocStruct,'costatename');
        statename=retrieveppdemodelinformation(ocStruct,'statename');
        statevector=cell2vectorstring([statename.value costatename.value]);
        specificdiscobjectivefunction=retrieveppdemodelinformation(ocStruct,'specificdiscobjectivefunction',arcidentifier,symkernel);
        data.type='mathcellchar';
        data.value=string2cell(ocmatjacobian(['[' specificdiscobjectivefunction.value ']'],statevector,symkernel),'vector');
        data.description='derivative of the discounted objective integrand with respect to the state and costate';

    case 'discobjectivefunctionDP'
        specificdiscobjectivefunction=retrieveppdemodelinformation(ocStruct,'specificdiscobjectivefunction',arcidentifier,symkernel);
        parametername=retrieveppdemodelinformation(ocStruct,'parametername',arcidentifier);
        variable=cell2vectorstring(parametername.value);
        data.type='mathcellchar';
        data.value=string2cell(ocmatjacobian(['[' specificdiscobjectivefunction.value ']'],variable,symkernel),'vector');
        data.description='derivative of the discounted objective integrand with respect to the state and costate';

    case 'hamiltonianfunction'
        % general Hamilton function including the "Lagrange terms" from
        % inequality contraints
        if isfield(ocStruct,'hamiltonianfunction') && isfield(ocStruct.hamiltonianfunction,'term')
            data.value=ocStruct.hamiltonianfunction.term;
        else
            statenum=retrieveppdemodelinformation(ocStruct,'statenum');
            costatename=retrieveppdemodelinformation(ocStruct,'costatename');
            % integrand of the objective value function without discount factor
            g=retrieveppdemodelinformation(ocStruct,'objectiveintegrand');
            % state dynamics
            dxdt=retrieveppdemodelinformation(ocStruct,'statedynamics');
            data.value=g.value;
            for ii=1:statenum.value
                data.value=[data.value '+' costatename.value{ii} '*(' dxdt.value{ii} ')'];
            end
        end
        data.type='mathchar';
        data.description='Hamilton function P=g+lamba*f';
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%
        % Control constraints
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
    case 'inequalitycontrolconstraintidentifier'
        if isfield(ocStruct.constraint,'function') && ...
                isfield(ocStruct.constraint.function,'control')
            data.value={ocStruct.constraint.function.control.identifier};
            data.type='cellchar';
        end
        
    case 'inequalitycontrolconstraint'
        % inequality constraints which explicitly include control variables
        % (and state variables).
        if isfield(ocStruct.constraint,'function') &&  ...
                isfield(ocStruct.constraint.function,'control')

            data.value={ocStruct.constraint.function.control.term};
            data.type='mathchar';
        end
        data.description='mixed inequality constraints';
        
    case 'inequalitycontrolconstraintnum'
        % number of inequality constraints explicitly including control
        % variables (and state variables).
        data.type='integer';
        if isfield(ocStruct.constraint,'function') && isfield(ocStruct.constraint.function,'control')
            data.value=length(ocStruct.constraint.function.control);
        else
            data.value=0;
        end
        data.description='number of mixed inequality constraints';


    case 'lagrangemultipliercontrolname'
        % variable name of the Lagrange multiplier(s) for inequality
        % (control) constraints as defined by the user in the
        % initialization file
        if isfield(ocStruct.variable,'lagrangemultcc')
            data.type='cellchar';
            data.value={ocStruct.variable.lagrangemultcc.name};
        end
        data.description='name of Lagrange multiplier';
        
    case 'constraintcombination'
        % allowed constraint combinations specified by the user in the
        % initialization file
        arcarg=arcidentifier2arcindex(arcidentifier);
        data.type='cellchar';
        data.value=regexp(strtrim(ocStruct.arc.constraintcombination{arcarg}),'\s*','split');
        data.description='allowed constraint combinations';

    case 'zerolmmc'
        % returns the index of the Lagrange multipliers corresponding to an
        % inactive constraint for a specific arc
        lagrangemultipliercontrolname=retrieveppdemodelinformation(ocStruct,'lagrangemultipliercontrolname');
        zerolmmcindex=retrieveppdemodelinformation(ocStruct,'zerolmmcindex',arcidentifier);
        data.value=lagrangemultipliercontrolname.value(zerolmmcindex.value);
        data.type='cellchar';
        data.description='name of the Lagrange multipliers corresponding to an inactive constraint';

    case 'zerolmmcindex'
        % returns the index of the Lagrange multipliers corresponding to an
        % inactive constraint for a specific arc
        inequalitycontrolconstraintidentifier=retrieveppdemodelinformation(ocStruct,'inequalitycontrolconstraintidentifier');
        controlconstraintnum=retrieveppdemodelinformation(ocStruct,'inequalitycontrolconstraintnum');
        constraintcombination=retrieveppdemodelinformation(ocStruct,'constraintcombination',arcidentifier);
        zerolmmcindex=[];
        if ~strcmp(constraintcombination.value,'[]')
            for jj=1:numel(inequalitycontrolconstraintidentifier.value)
                testexpr=regexp(constraintcombination.value,['\<' inequalitycontrolconstraintidentifier.value{jj} '\>'],'ONCE');
                if isempty([testexpr{:}])
                    zerolmmcindex=[zerolmmcindex jj];
                end
            end
        else
            zerolmmcindex=1:controlconstraintnum.value;
        end
        data.value=zerolmmcindex;
        data.type='integer';
        data.description='index of the Lagrange multipliers corresponding to an inactive constraint';

    case 'nonzerolmmc'
        % returns the name of the Lagrange multipliers corresponding to an
        % active constraint for a specific arc
        lagrangemultipliercontrolname=retrieveppdemodelinformation(ocStruct,'lagrangemultipliercontrolname');
        nonzerolmmcindex=retrieveppdemodelinformation(ocStruct,'nonzerolmmcindex',arcidentifier);
        data.value=lagrangemultipliercontrolname.value(nonzerolmmcindex.value);
        data.type='cellchar';
        data.description='name of the Lagrange multipliers corresponding to an active constraint';

    case 'nonzerolmmcindex'
        % returns the index of the Lagrange multipliers corresponding to an
        % active constraint for a specific arc
        inequalitycontrolconstraintnum=retrieveppdemodelinformation(ocStruct,'inequalitycontrolconstraintnum');
        zerolmmcindex=retrieveppdemodelinformation(ocStruct,'zerolmmcindex',arcidentifier);
        data.value=setdiff(1:inequalitycontrolconstraintnum.value,zerolmmcindex.value);
        data.type='integer';
        data.description='index of the Lagrange multipliers corresponding to an active constraint';

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Maximizing Condition
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
    case 'optimizationtype'
        data.type='char';
        data.value=ocStruct.optimizationtype;
        data.description='allowed constraint combinations';

    case 'pontryaginfunction'
        % general Hamilton function including the "Lagrange terms" from
        % inequality contraints
        if isfield(ocStruct,'pontryaginfunction') && isfield(ocStruct.pontryaginfunction,'term')
            data.value=ocStruct.pontryaginfunction.term;
        else
            statenum=retrieveppdemodelinformation(ocStruct,'statenum');
            costatename=retrieveppdemodelinformation(ocStruct,'costatename');
            inequalitycontrolconstraint=retrieveppdemodelinformation(ocStruct,'inequalitycontrolconstraint');
            inequalitycontrolconstraintnum=retrieveppdemodelinformation(ocStruct,'inequalitycontrolconstraintnum');
            lmmc=retrieveppdemodelinformation(ocStruct,'lagrangemultipliercontrolname');
            % integrand of the objective value function without discount factor
            g=retrieveppdemodelinformation(ocStruct,'objectiveintegrand');
            % state dynamics
            dxdt=retrieveppdemodelinformation(ocStruct,'statedynamics');
            data.value=g.value;
            for ii=1:statenum.value
                data.value=[data.value '+' costatename.value{ii} '*(' dxdt.value{ii} ')'];
            end
            for ii=1:inequalitycontrolconstraintnum.value
                data.value=[data.value '+' lmmc.value{ii} '*(' inequalitycontrolconstraint.value{ii} ')'];
            end
        end
        data.type='mathchar';
        data.description='general Hamilton function P=H+mu*mc';
        

    case 'pontryaginfunctionDx'
        % derivative of the Pontryaginfunction with respect to the state
        % variable(s)
        if isfield(ocStruct,'pontryaginfunction') && isfield(ocStruct.pontryaginfunction,'derivative') ...
                && isfield(ocStruct.pontryaginfunction.derivative,'Dx')
            data.value=ocStruct.pontryaginfunction.derivative.Dx.term;
        else
            statename=retrieveppdemodelinformation(ocStruct,'statename');
            statevector=cell2vectorstring(statename.value);
            pontryaginfunction=retrieveppdemodelinformation(ocStruct,'pontryaginfunction');
            data.value=string2cell(ocmatjacobian(['[' pontryaginfunction.value ']'],statevector,symkernel),'vector');
        end
        data.type='mathcellchar';
        data.description='derivative of the Pontryaginfunction with respect to the state';

    case 'pontryaginfunctionDX'
        % derivative of the Pontryaginfunction with respect to the state
        % AND costate variable(s)
        costatename=retrieveppdemodelinformation(ocStruct,'costatename');
        statename=retrieveppdemodelinformation(ocStruct,'statename');
        statevector=cell2vectorstring([statename.value costatename.value]);
        pontryaginfunction=retrieveppdemodelinformation(ocStruct,'specificpontryaginfunction',arcidentifier,symkernel);
        data.value=string2cell(ocmatjacobian(['[' pontryaginfunction.value ']'],statevector,symkernel),'vector');
        data.type='mathcellchar';
        data.description='derivative of the Pontryaginfunction with respect to the state and costate';


    case 'pontryaginfunctionDu'
        % derivative of the Pontryaginfunction with respect to the control
        % variable(s)
        if isfield(ocStruct,'pontryaginfunction') && isfield(ocStruct.pontryaginfunction,'derivative') ...
                && isfield(ocStruct.pontryaginfunction.derivative,'Du')
            data.value=ocStruct.pontryaginfunction.derivative.Du.term;
        else
            controlname=retrieveppdemodelinformation(ocStruct,'controlname');
            controlvector=cell2vectorstring(controlname.value);
            pontryaginfunction=retrieveppdemodelinformation(ocStruct,'pontryaginfunction');
            data.value=string2cell(ocmatjacobian(['[' pontryaginfunction.value ']'],controlvector,symkernel),'vector');
        end
        data.type='mathcellchar';
        data.description='derivative of the Pontryaginfunction with respect to the control';

    case 'pontryaginfunctionDx2'
        % second order derivative of the Pontryaginfunction with respect to
        % the state variable(s)
        if isfield(ocStruct,'pontryaginfunction') && isfield(ocStruct.pontryaginfunction,'derivative') ...
                && isfield(ocStruct.pontryaginfunction.derivative,'Dx2')
            data.value=ocStruct.pontryaginfunction.derivative.Dx2.term;
        else
            statename=retrieveppdemodelinformation(ocStruct,'statename');
            statevector=cell2vectorstring(statename.value);
            pontryaginfunction=retrieveppdemodelinformation(ocStruct,'pontryaginfunction');
            data.value=string2cell(ocmathessian(pontryaginfunction.value,statevector,symkernel),'matrix');
        end
        data.type='mathcellchar';
        data.description='second order derivative of the Pontryaginfunction with respect to the state';

    case 'pontryaginfunctionDu2'
        % second order derivative of the Pontryaginfunction with respect to
        % the control variable(s)
        if isfield(ocStruct,'pontryaginfunction') && isfield(ocStruct.pontryaginfunction,'derivative') ...
                && isfield(ocStruct.pontryaginfunction.derivative,'Du2')
            data.value=ocStruct.pontryaginfunction.derivative.Du2.term;
        else
            controlname=retrieveppdemodelinformation(ocStruct,'controlname');
            controlvector=cell2vectorstring(controlname.value);
            pontryaginfunction=retrieveppdemodelinformation(ocStruct,'pontryaginfunction');
            data.value=string2cell(ocmathessian(pontryaginfunction.value,controlvector,symkernel),'matrix');
        end
        data.type='mathcellchar';
        data.description='second order derivative of the Pontryaginfunction with respect to the control';

    case {'pontryaginfunctionDuDx','pontryaginfunctionDxDu'}
        % second order derivative of the Pontryaginfunction with respect to
        % the control variable(s)
        if isfield(ocStruct,'pontryaginfunction') && isfield(ocStruct.pontryaginfunction,'derivative') ...
                && isfield(ocStruct.pontryaginfunction.derivative,'DuDx')
            data.value=ocStruct.pontryaginfunction.derivative.DuDx.term;
        else
            controlname=retrieveppdemodelinformation(ocStruct,'controlname');
            controlvector=cell2vectorstring(controlname.value);
            statename=retrieveppdemodelinformation(ocStruct,'statename');
            statevector=cell2vectorstring(statename.value);
            pontryaginfunction=retrieveppdemodelinformation(ocStruct,'pontryaginfunction');
            J=removematrixstring(ocmatjacobian(['[' pontryaginfunction.value ']'],controlvector,symkernel));
            data.value=string2cell(ocmatjacobian(J,statevector,symkernel),'matrix');
        end
        data.type='mathcellchar';
        data.description='second order derivative of the Pontryaginfunction with respect to the control';
        
    case 'pontryaginfunction4arc'
        % returns the specific Pontryagin functions for every considered
        % constraint combination. It is assumed that the constraints are
        % given by ODEs (statedynamics) and inequality constraints
        % (control-state constraints)
        pontryaginfunction=retrieveppdemodelinformation(ocStruct,'pontryaginfunction');
        zerolmmcindex=retrieveppdemodelinformation(ocStruct,'zerolmmcindex',arcidentifier);
        lagrangemultipliercontrolname=retrieveppdemodelinformation(ocStruct,'lagrangemultipliercontrolname');
        lmmcequalzero='';
        ii=0;
        for idx=zerolmmcindex.value
            ii=ii+1;
            lmmcequalzero{ii}=[lagrangemultipliercontrolname.value{idx} '=0'];
        end
        zerolmmcvector=cell2vectorstring(lmmcequalzero);
        if ~isempty(zerolmmcvector)
            pontryaginfunction.value=ocmatsubs(pontryaginfunction.value,zerolmmcvector,symkernel);
        end
        data.value=pontryaginfunction.value;
        data.type='mathchar';
        data.description='specific Hamilton function for actual constraint combination';

    case 'arcwitheequalpontryaginfunction'
        % in cases where the Hamiltonian maximizing condition yields
        % multiple result different solutions are identified by different
        % arc identifiers. 'arcwitheequalpontryaginfunction' returns the
        % arc idnetifiers which correspond to the same specification of
        % the Pontryagin function (generalized Hamilton function including
        % terms for (inequality) constraints)
        uniqueconstraintcombination=unique(ocStruct.arc.constraintcombination);
        arcidentifier=ocStruct.arc.identifier;
        if numel(uniqueconstraintcombination)==numel(ocStruct.arc.constraintcombination)
            for ii=1:numel(uniqueconstraintcombination)
                data.value{ii}=arcidentifier(ii);
            end
        else
            for ii=1:numel(uniqueconstraintcombination)
                data.value{ii}=arcidentifier(strcmp(ocStruct.arc.constraintcombination,uniqueconstraintcombination{ii}));
            end
        end
        data.type='cellchar';
        data.description='arcidentifer corresponding to same Pontryaginfunction';

    case 'maximizingvariable'
        % dependent variable names of the Hamiltonian maximizing condition
        controlname=retrieveppdemodelinformation(ocStruct,'controlname');
        lagrangemultipliercontrolname=retrieveppdemodelinformation(ocStruct,'lagrangemultipliercontrolname');
        data.value=[controlname.value lagrangemultipliercontrolname.value];
        data.type='cellchar';
        data.description='variable names of the Hamiltonian maximizing condition';

    case 'solutionindex'
        % if more than one solution of Hamiltonian maximizing condition
        % exists this index determines which of the solutions skould be kept
        % by default. During the initialization porcess the user is
        % explicitly asked to confirm this selection
        uniqueconstraintcombination=unique(ocStruct.arc.constraintcombination);
        controlvaluecombination=ocStruct.arc.controlvaluecombination;
         if numel(uniqueconstraintcombination)==numel(ocStruct.arc.constraintcombination)
            data.value=controlvaluecombination;
        else
            for ii=1:numel(uniqueconstraintcombination)
                data.value{ii}=controlvaluecombination(strcmp(ocStruct.arc.constraintcombination,uniqueconstraintcombination{ii}));
            end
        end
        data.description='used index for multiple Hamilton maximizing condition solutions';

    case 'explicitnonlinearcontrol'
        % explicit nonlinear control variables for a specific arc
        struct=[ocStruct.variable.control.(arcidentifier2field(arcidentifier))];
        for ii=1:length(struct)
            data.value{ii}=ocStruct.variable.control(ii).name([struct(ii).property.explicit]==1 & [struct(ii).property.linear]==0);
        end
        data.type='cellchar';
        data.description='explicit nonlinear control variable name';

    case 'explicitnonlinearcontrolindex'
        % explicit nonlinear control variables for a specific arc
        struct=[ocStruct.variable.control.(arcidentifier2field(arcidentifier))];
        for ii=1:length(struct)
            data.value(ii)=([struct(ii).property.explicit]==1 & [struct(ii).property.linear]==0);
        end
        data.type='integer vector';
        data.description='explicit nonlinear control variable name';

    case 'maximizingexplicitvariable'
        % variable names which are used to derive the first order necessary
        % conditions of the Hamiltonian maximizing condition and are
        % calculated explicitly
        arc=arcidentifier2field(arcidentifier);
        if isfield(ocStruct,'foc') &&isfield(ocStruct.foc,'generalinformation') && ...
                isfield(ocStruct.foc.generalinformation,arc) && ...
                isfield(ocStruct.foc.generalinformation.(arc),'maximizingexplicitvariable')
            data.value=ocStruct.foc.generalinformation.(arc).maximizingexplicitvariable.name;
        else
            explicitnonlinearcontrol=retrieveppdemodelinformation(ocStruct,'explicitnonlinearcontrol',arcidentifier);
            nonzerolmmc=retrieveppdemodelinformation(ocStruct,'nonzerolmmc',arcidentifier);
            data.value=[explicitnonlinearcontrol.value nonzerolmmc.value];
        end
        data.type='cellchar';
        data.description='explicit variable names for the first order necessary conditions of the Hamiltonian maximizing condition';
        
    case 'maximizingregularnonlinearsolution'
        % Pontryagin's maximumprinciple is used to derive explicit formulas
        % of control and Lagrange multiplier values satisfying the first
        % order necessary conditions
        solution=[];
        totalmaximizingvariable=retrieveppdemodelinformation(ocStruct,'maximizingvariable');
        maximizingexplicitvariable=retrieveppdemodelinformation(ocStruct,'maximizingexplicitvariable',arcidentifier);
        zerolmmc=retrieveppdemodelinformation(ocStruct,'zerolmmc',arcidentifier);
        %nonzerolmmc=retrieveppdemodelinformation(ocStruct,'nonzerolmmc',arcidentifier);
        controlname=retrieveppdemodelinformation(ocStruct,'controlname');
        lagrangemultipliercontrolname=retrieveppdemodelinformation(ocStruct,'lagrangemultipliercontrolname');
        controlnum=retrieveppdemodelinformation(ocStruct,'controlnum');
        inequalitycontrolconstraintnum=retrieveppdemodelinformation(ocStruct,'inequalitycontrolconstraintnum');
        pontryaginfunction4arc=retrieveppdemodelinformation(ocStruct,'pontryaginfunction4arc',arcidentifier,symkernel);
        controldependence=retrieveppdemodelinformation(ocStruct,'controldependence');
        spacearg=retrieveppdemodelinformation(ocStruct,'space');
        controltimedependence=strcmp(controldependence.value,'0');

        if ~isempty(maximizingexplicitvariable.value)
            equationvariable=cell2vectorstring(maximizingexplicitvariable.value);
            equation=removematrixstring(ocmatjacobian(['[' pontryaginfunction4arc.value ']'],equationvariable,symkernel));
            % ocmatsolve is an adaptation of the native MATLAB command
            % solve for the symbolic toolbox relying
            if ~strcmp(equation,'[]')
                solution=ocmatsolve(equation,equationvariable,symkernel);
                if isempty(solution)
                    %ocmatmsg('Explicit solution could not be found.\n')
                    return
                end
            end
            if ~isempty(solution)
                for jj=1:numel(solution)
                    for ii=1:numel(zerolmmc.value)
                        % the Lagrange multipliers being zero (inactive
                        % constraint) are added
                        solution(jj).(zerolmmc.value{ii})='0';
                    end
                end
            else
                for ii=1:numel(zerolmmc.value)
                    % the Lagrange multipliers being zero (inactive
                    % constraint) are added
                    solution.(zerolmmc.value{ii})='0';
                end
            end
        else
            solution=[];
        end
        solution=orderfields(solution,totalmaximizingvariable.value);
        for ii=1:numel(solution)
            for jj=1:controlnum.value
                if controltimedependence(jj)
                    solution(ii).control{jj}=['int_' spacearg.value '(' solution(ii).(controlname.value{jj}) ')'];
                else
                     solution(ii).control{jj}=solution(ii).(controlname.value{jj});
               end
            end
            for jj=1:inequalitycontrolconstraintnum.value
                solution(ii).lagrangemultcc{jj}=solution(ii).(lagrangemultipliercontrolname.value{jj});
            end
        end
        data.value=solution;
        data.type='struct';
        data.description='solution(s) derived from the foc of the Hamiltonian maximizing condition';

    case 'optimalvalue'
        % the solutions from the Hamiltonian maximizing condition have
        % already to be stored in the ocStruct structure.
        controlname=retrieveppdemodelinformation(ocStruct,'controlname');
        controlnum=retrieveppdemodelinformation(ocStruct,'controlnum');
        lagrangemultipliercontrolname=retrieveppdemodelinformation(ocStruct,'lagrangemultipliercontrolname');
        inequalitycontrolconstraintnum=retrieveppdemodelinformation(ocStruct,'inequalitycontrolconstraintnum');

        for ii=1:controlnum.value
            value.(controlname.value{ii})=ocStruct.foc.value.control.(arcidentifier2field(arcidentifier)).term{ii};
        end
        for ii=1:inequalitycontrolconstraintnum.value
            value.(lagrangemultipliercontrolname.value{ii})=ocStruct.foc.value.lagrangemultcc.(arcidentifier2field(arcidentifier)).term{ii};
        end
        data.type='structmathchar';
        data.value=value;
        data.description='optimal control values as derived from the Hamiltonian maximizing condition';
        
        %%%%%%%%%%%%%%%%%%%%%%%%
        % Adjoint System
        %%%%%%%%%%%%%%%%%%%%%%%%

    case 'canonicalsystem'
        % the canonical system consisting of the state dynamics, adjoint
        % system and possible algebraic equations for implicitly given
        % control values.
        adjointsystem=retrieveppdemodelinformation(ocStruct,'adjointsystem','',symkernel);
        statedynamics=retrieveppdemodelinformation(ocStruct,'statedynamics');
        dxdt=[statedynamics.value adjointsystem.value];
        data.value=dxdt;
        data.type='mathcellchar';
        data.description='state dynamics, adjoint system';

    case 'adjointsystem'
        % the adjoint system, where each cell contains one costate dynamics
        arc=arcidentifier2field(arcidentifier);
        if 0%isfield(ocStruct,'foc') && isfield(ocStruct.foc,'adjointsystem') ...
                %&& isfield(ocStruct.foc.adjointsystem,'dynamics') ...
                %&& isfield(ocStruct.foc.adjointsystem.dynamics,arc) ...
            data.value=ocStruct.foc.adjointsystem.dynamics.(arc).term;
        else
            pontryaginfunctionDx=retrieveppdemodelinformation(ocStruct,'pontryaginfunctionDx','',symkernel);
            statenum=retrieveppdemodelinformation(ocStruct,'statenum');
            costatename=retrieveppdemodelinformation(ocStruct,'costatename');
            discountratevariable=retrieveppdemodelinformation(ocStruct,'discountratevariable');
            for ii=1:statenum.value
                data.value{ii}=[discountratevariable.value '*' costatename.value{ii} '-(' pontryaginfunctionDx.value{ii} ')'];
            end
        end
        data.type='mathcellchar';
        data.description='costate dynamics';

    case 'adjointsystemboundarycondition'
        % the adjoint system, where each cell contains one costate dynamics
        arc=arcidentifier2field(arcidentifier);
        statenum=retrieveppdemodelinformation(ocStruct,'statenum');
        if 0
        else
            stateboundarycondition=retrieveppdemodelinformation(ocStruct,'stateboundarycondition');
            for ii=1:statenum.value
                switch stateboundarycondition.value{ii}.type
                    case 'Neumann1D'
                        stateboundarycondition.value{ii}.term='0';
                        data.value{ii}=stateboundarycondition.value{ii};
                end
            end
            
        end
        data.type='cell of structs';
        data.description='spatial boundary conditions for the costate(s)';

    case 'adjointsystemconvection'
        % the adjoint system, where each cell contains one costate dynamics
        arc=arcidentifier2field(arcidentifier);
        if isfield(ocStruct,'foc') && isfield(ocStruct.foc,'adjointsystem') ...
                && isfield(ocStruct.foc.adjointsystem,'dynamics') ...
                && isfield(ocStruct.foc.adjointsystem.dynamics,arc) ...
                data.value=ocStruct.foc.adjointsystem.dynamics.(arc).diffusion;
        else
            statedynamicsconvection=retrieveppdemodelinformation(ocStruct,'statedynamicsconvection');
            statedynamicsconvectionvariable=retrieveppdemodelinformation(ocStruct,'statedynamicsconvectionvariable');
            statenum=retrieveppdemodelinformation(ocStruct,'statenum');
            statename=retrieveppdemodelinformation(ocStruct,'statename');
            costatename=retrieveppdemodelinformation(ocStruct,'costatename');
            for ii=1:statenum.value
                if ~isempty(statedynamicsconvectionvariable.value)
                    fidx=find(strcmp(statename.value,statedynamicsconvectionvariable.value{ii}));
                    costatedynamicsconvectionvariable=costatename.value{fidx};
                    if strcmp(statedynamicsconvection.value{ii}(1),'-')
                        costatedynamicsconvection=statedynamicsconvection.value{ii}(2:end);
                    else
                        costatedynamicsconvection=['-(' statedynamicsconvection.value{ii} ')'];
                    end
                    data.value.variable{ii}=costatedynamicsconvectionvariable;
                    data.value.term{ii}=costatedynamicsconvection;
                else
                    data.value=[];
                end
            end
        end
        data.type='mathcellchar';
        data.description='costate dynamics';

    case 'adjointsystemdiffusion'
        % the adjoint system, where each cell contains one costate dynamics
        arc=arcidentifier2field(arcidentifier);
        if 0 %isfield(ocStruct,'foc') && isfield(ocStruct.foc,'adjointsystem') ...
              %  && isfield(ocStruct.foc.adjointsystem,'dynamics') ...
              %  && isfield(ocStruct.foc.adjointsystem.dynamics,arc) ...
                data.value=ocStruct.foc.adjointsystem.dynamics.(arc).diffusion;
        else
            statedynamicsdiffusion=retrieveppdemodelinformation(ocStruct,'statedynamicsdiffusion');
            statedynamicsdiffusionvariable=retrieveppdemodelinformation(ocStruct,'statedynamicsdiffusionvariable');
            statenum=retrieveppdemodelinformation(ocStruct,'statenum');
            statename=retrieveppdemodelinformation(ocStruct,'statename');
            costatename=retrieveppdemodelinformation(ocStruct,'costatename');
            for ii=1:statenum.value
                if ~isempty(statedynamicsdiffusionvariable.value)
                    fidx=find(strcmp(statename.value,statedynamicsdiffusionvariable.value{ii}));
                    costatedynamicsdiffusionvariable=costatename.value{fidx};
                    if strcmp(statedynamicsdiffusion.value{ii}(1),'-')
                        costatedynamicsdiffusion=statedynamicsdiffusion.value{ii}(2:end);
                    else
                        costatedynamicsdiffusion=['-(' statedynamicsdiffusion.value{ii} ')'];
                    end
                    data.value{ii}.variable=costatedynamicsdiffusionvariable;
                    data.value{ii}.term=costatedynamicsdiffusion;
                else
                    data.value=[];
                end
            end
        end
        data.type='mathcellchar';
        data.description='costate dynamics';
        
    case 'specificstatedynamics'
        % the canonical system consisting of the state dynamics, adjoint
        % system and possible algebraic equations for implicitly given
        % control values.
        statedynamics=retrieveppdemodelinformation(ocStruct,'statedynamics');
        optimalvalue=retrieveppdemodelinformation(ocStruct,'optimalvalue',arcidentifier);
        maximizingvariable=retrieveppdemodelinformation(ocStruct,'maximizingvariable');
        dxdt=statedynamics.value;

        for ii=1:numel(maximizingvariable.value)
            for jj=1:numel(dxdt)
                dxdt{jj}=ocmatsubs(dxdt{jj},[maximizingvariable.value{ii} '=' optimalvalue.value.(maximizingvariable.value{ii})],symkernel);
            end
        end
        data.value=dxdt;
        data.type='mathcellchar';
        data.description='state dynamics';

    case 'specificadjointsystem'
        % the canonical system consisting of the state dynamics, adjoint
        % system and possible algebraic equations for implicitly given
        % control values.
        adjointsystem=retrieveppdemodelinformation(ocStruct,'adjointsystem','',symkernel);
        optimalvalue=retrieveppdemodelinformation(ocStruct,'optimalvalue',arcidentifier);
        maximizingvariable=retrieveppdemodelinformation(ocStruct,'maximizingvariable');
        dldt=adjointsystem.value;

        for ii=1:numel(maximizingvariable.value)
            for jj=1:numel(dldt)
                dldt{jj}=ocmatsubs(dldt{jj},[maximizingvariable.value{ii} '=' optimalvalue.value.(maximizingvariable.value{ii})],symkernel);
            end
        end
        data.value=dldt;
        data.type='mathcellchar';
        data.description='costate dynamics';

    case 'specificcanonicalsystem'
        % the canonical system consisting of the state dynamics, adjoint
        % system and possible algebraic equations for implicitly given
        % control values.
        adjointsystem=retrieveppdemodelinformation(ocStruct,'adjointsystem','',symkernel);
        statedynamics=retrieveppdemodelinformation(ocStruct,'statedynamics');
        optimalvalue=retrieveppdemodelinformation(ocStruct,'optimalvalue',arcidentifier);
        maximizingvariable=retrieveppdemodelinformation(ocStruct,'maximizingvariable');
        dxdt=[statedynamics.value adjointsystem.value];
        for ii=1:numel(maximizingvariable.value)
            for jj=1:numel(dxdt)
                if ~strcmp(dxdt{jj},'0')%~(dxdt{jj}==sym('0'))
                    dxdt{jj}=ocmatsubs(dxdt{jj},[maximizingvariable.value{ii} '=' optimalvalue.value.(maximizingvariable.value{ii})],symkernel);
                end
            end
        end
        data.value=dxdt;
        data.type='mathcellchar';
        data.description='state dynamics, adjoint system (and algebraic equations for implicit control values)';

    case 'canonicalsystemjacobian'
        % Jacobian of the canonical system with respect to state, costate
        % and implicit control/Lagrange multiplier variables
        data.value=ocStruct.foc.canonicalsystem.derivative.(arcidentifier2field(arcidentifier)).DX.term;
        data.type='mathcellchar';
        data.description='Jacobian of the canonical system with respect to state, costate and implicit control/Lagrange multiplier variables';
        
        %%%%%%%%%%%%%%%%%%%%%%%%%
        % Model Parameter
        %%%%%%%%%%%%%%%%%%%%%%%%%
    case 'parametername'
        % variable name(s) of the parameter(s)
        data.type='cellchar';
        data.value=fieldnames(ocStruct.parameter.variable);
        data.description='name of parameter variable';

    case 'parametervalue'
        % vector of user provided parameter values
        data.type='double';
        parametername=retrieveodemodelinformation(ocStruct,'parametername');
        data.value=zeros(1,numel(parametername.value));
        for ii=1:numel(parametername.value)
            data.value(ii)=ocStruct.parameter.variable.(parametername.value{ii});
        end
        data.description='vector of user provided parameter values';
        
    case 'discountfactor'
        % discount variable used in the discount factor of the objective
        % value function
        if isfield(ocStruct.objective,'integral')
            data.type='mathchar';
            data.value=ocStruct.objective.integral.discountfactor.term;
            data.description='discount factor of the objective value function';
        else
            data.type='mathchar';
            data.value=[];
            data.description='discount factor of the objective value function';
        end
        data.type='mathchar';
        data.description='discount factor of the objective value function';
    case 'discountrate'
        data.type='double';
        data.value=ocStruct.objective.integral.discountrate;
        data.description='value of the discountrate';
        
    case 'algebraicequationnum'
        data.type='integer';
        data.value=0;
        data.description='number of algebraic equations';

    case 'variablename'
        % general variable name(s) of the optimal control problem
        data.type='cellchar';
        data.value=fieldnames(ocStruct.variable);
        data.description='general variable name';

    case 'equationvariablename'
        % dependent variables of the canonical system for a specific arc
        data.type='cellchar';
        statename=retrieveppdemodelinformation(ocStruct,'statename');
        data.value=statename.value;
        data.description='dependent variables of the canonical system';

    case 'statenum'
        % the number of state(s)
        data.type='integer';
        data.value=length(ocStruct.variable.state);
        data.description='number of state';


    case 'controlnum'
        % the number of state(s)
        data.type='integer';
        data.value=length(ocStruct.variable.control);
        data.description='number of controls';

    case 'arcnum'
        % the number of different arcs, corresponding to different
        % functional specifications of the canonical system
        data.type='integer';
        data.value=ocStruct.arc.num;
        data.description='number of different arcs';

    case 'odedim'
        % the number of state(s)
        data.type='integer';
        data.value=ocStruct.variable.state.num;
        data.description='number of ODEs for the canonical system';

    case 'dynamics'
        % ODEs describing the evolution of the state(s)
        data.type='mathchar';
        data.value=ocStruct.constraint.derivative.state.term;
        data.description='ODEs of the dynamics';

    case 'dynamicsjacobian'
        % Jacobian of the canonical system with respect to state, costate
        % and implicit control/Lagrange multiplier variables
        dynamics=retrieveodemodelinformation(ocStruct,'dynamics');
        dxdt=cell2vectorstring(dynamics.value);
        statename=retrieveodemodelinformation(ocStruct,'statename');
        variable=cell2vectorstring(statename.value);
        data.value=string2cell(removematrixstring(ocmatjacobian(dxdt,variable,symkernel)),'matrix');
        data.type='mathcellchar';
        data.description='Jacobian of the dynamics with respect to state';

    case 'dynamicsparameterjacobian'
        % Jacobian of the canonical system with respect to state, costate
        % and implicit control/Lagrange multiplier variables
        dynamics=retrieveodemodelinformation(ocStruct,'dynamics');
        dxdt=cell2vectorstring(dynamics.value);
        parametername=retrieveodemodelinformation(ocStruct,'parametername');
        variable=cell2vectorstring(parametername.value);
        data.value=string2cell(removematrixstring(ocmatjacobian(dxdt,variable,symkernel)),'matrix');
        data.type='mathcellchar';
        data.description='Jacobian of the dynamics with respect to parametervalues';
end



function cellstring=string2cell(string,matrixtype)
% each cell for matrixtype 'matrix' consists of a row of the matrix
% each cell for matrixtype 'vector' consists of an entry of the vector

cellstring=[];
string=removematrixstring(string);
if ~strcmp(string([1 end]),'[]')
    cellstring{1}=string;
    return
end
switch matrixtype
    case 'vector'
        cellstring=regexp(string(2:end-1),',','split');
    case 'matrix'
        cellstring=regexp(string(2:end-1),'],(\ )*[','split');
    case 'charmatrix'
        cellstring=regexp(string(2:end-1),'],[','split');
end

function symval=mycell2sym(cellvec,matrixtype)
% each cell for matrixtype 'matrix' consists of a row of the matrix
% each cell for matrixtype 'vector' consists of an entry of the vector

symval=sym([]);
l=length(cellvec);
switch matrixtype
    case 'vector'
        stringval=['[' deblank(cellvec{1})];
        for ii=2:l
            stringval=[stringval ';' deblank(cellvec{ii})];
        end
        symval=sym([stringval ']']);
    case 'matrix'
        stringval=['[' deblank(cellvec{1})];
        for ii=2:l
            stringval=[stringval ';' deblank(cellvec{ii})];
        end
        symval=sym([stringval ']']);
end

