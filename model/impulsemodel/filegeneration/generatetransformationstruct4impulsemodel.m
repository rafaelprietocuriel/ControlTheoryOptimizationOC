function transStruct=generatetransformationstruct4impulsemodel(ocStruct,transStruct,varname,overwrite,varargin)
%
% in this structure the information for the transformation of optimal
% control specific variables to default names, which are used in the
% automatically generated files are stored.
% 
if nargin<=3
    overwrite=0;
end

if ~overwrite
    if isfield(transStruct,varname)
        return
    end
end

switch varname
    case 'independent'
        transStruct.(varname).arcdependent=0;
        transStruct.(varname).name=getbasicname(varname);
    case 'state'
        basename=getbasicname('dependent');
        statenum=retrieveimpulsemodelinformation(ocStruct,'statenum');
        transStruct.(varname).arcdependent=0;
        transStruct.(varname).basicname=basename;
        for ii=1:statenum.value
            transStruct.(varname).name{ii}=[basename int2str(ii)];
            transStruct.(varname).arrayname{ii}=[basename '(' int2str(ii) ')'];
            transStruct.(varname).vectorizedname{ii}=[basename '(' int2str(ii) ',:)'];
        end
    case 'statel'
        basename=[getbasicname('dependent') 'L'];
        statenum=retrieveimpulsemodelinformation(ocStruct,'statenum');
        transStruct.(varname).arcdependent=0;
        transStruct.(varname).basicname=basename;
        for ii=1:statenum.value
            transStruct.(varname).name{ii}=[basename int2str(ii)];
            transStruct.(varname).arrayname{ii}=[basename '(' int2str(ii) ')'];
            transStruct.(varname).vectorizedname{ii}=[basename '(' int2str(ii) ',:)'];
        end
    case 'stater'
        basename=[getbasicname('dependent') 'R'];
        statenum=retrieveimpulsemodelinformation(ocStruct,'statenum');
        transStruct.(varname).arcdependent=0;
        transStruct.(varname).basicname=basename;
        for ii=1:statenum.value
            transStruct.(varname).name{ii}=[basename int2str(ii)];
            transStruct.(varname).arrayname{ii}=[basename '(' int2str(ii) ')'];
            transStruct.(varname).vectorizedname{ii}=[basename '(' int2str(ii) ',:)'];
        end
    case 'costate'
        basename=getbasicname('dependent');
        costatenum=retrieveimpulsemodelinformation(ocStruct,'costatenum');
        transStruct.(varname).arcdependent=0;
        transStruct.(varname).basicname=basename;
        for ii=1:costatenum.value
            costatecoord=ii+ocStruct.variable.state.num;
            transStruct.(varname).name{ii}=[basename int2str(costatecoord)];
            transStruct.(varname).arrayname{ii}=[basename '(' int2str(costatecoord) ')'];
            transStruct.(varname).vectorizedname{ii}=[basename '(' int2str(costatecoord) ',:)'];
        end
    case 'costatel'
        basename=[getbasicname('dependent') 'L']; 
        costatenum=retrieveimpulsemodelinformation(ocStruct,'costatenum');
        transStruct.(varname).arcdependent=0;
        transStruct.(varname).basicname=basename;
        for ii=1:costatenum.value
            costatecoord=ii+ocStruct.variable.state.num;
            transStruct.(varname).name{ii}=[basename int2str(costatecoord)];
            transStruct.(varname).arrayname{ii}=[basename '(' int2str(costatecoord) ')'];
            transStruct.(varname).vectorizedname{ii}=[basename '(' int2str(costatecoord) ',:)'];
        end
    case 'costater'
        basename=[getbasicname('dependent') 'R']; 
        costatenum=retrieveimpulsemodelinformation(ocStruct,'costatenum');
        transStruct.(varname).arcdependent=0;
        transStruct.(varname).basicname=basename;
        for ii=1:costatenum.value
            costatecoord=ii+ocStruct.variable.state.num;
            transStruct.(varname).name{ii}=[basename int2str(costatecoord)];
            transStruct.(varname).arrayname{ii}=[basename '(' int2str(costatecoord) ')'];
            transStruct.(varname).vectorizedname{ii}=[basename '(' int2str(costatecoord) ',:)'];
        end
    case 'control'
        basename=getbasicname(varname);
        controlnum=retrieveimpulsemodelinformation(ocStruct,'controlnum');
        arcnum=retrieveimpulsemodelinformation(ocStruct,'arcnum');
        arcid=retrieveimpulsemodelinformation(ocStruct,'arcidentifier');
        transStruct.(varname).arcdependent=1;
        transStruct.(varname).basicname=basename;
        if controlnum.value
            for ii=1:controlnum.value
                transStruct.(varname).name{ii}=[basename int2str(ii)];
                transStruct.(varname).arrayname{ii}=[basename '(' int2str(ii) ')'];
                transStruct.(varname).vectorizedname{ii}=[basename '(' int2str(ii) ',:)'];
            end
            for ii=1:arcnum.value
                arc=arcidentifier2field(arcid.value{ii});
                for jj=1:controlnum.value
                    transStruct.(varname).(arc).name{jj}=[basename int2str(jj)];
                    transStruct.(varname).(arc).arrayname{jj}=[basename '(' int2str(jj) ')'];
                    transStruct.(varname).(arc).vectorizedname{jj}=[basename '(' int2str(jj) ',:)'];
                end
            end
        else
            transStruct.(varname).name=[];
            transStruct.(varname).arrayname=[];
            transStruct.(varname).vectorizedname=[];
            for ii=1:arcnum.value
                arc=arcidentifier2field(arcid.value{ii});
                transStruct.(varname).(arc).name=[];
                transStruct.(varname).(arc).arrayname=[];
                transStruct.(varname).(arc).vectorizedname=[];
            end
        end

    case 'impulsecontrol'
        basename=getbasicname(varname);
        controlnum=retrieveimpulsemodelinformation(ocStruct,'icontrolnum');
        arcnum=retrieveimpulsemodelinformation(ocStruct,'arcnum');
        arcid=retrieveimpulsemodelinformation(ocStruct,'arcidentifier');
        transStruct.(varname).arcdependent=1;
        transStruct.(varname).basicname=basename;
        for ii=1:controlnum.value
            transStruct.(varname).name{ii}=[basename int2str(ii)];
            transStruct.(varname).arrayname{ii}=[basename '(' int2str(ii) ')'];
            transStruct.(varname).vectorizedname{ii}=[basename '(' int2str(ii) ',:)'];
        end
        for ii=1:arcnum.value
            arc=arcidentifier2field(arcid.value{ii});
            for jj=1:controlnum.value
                transStruct.(varname).(arc).name{jj}=[basename int2str(jj)];
                transStruct.(varname).(arc).arrayname{jj}=[basename '(' int2str(jj) ')'];
                transStruct.(varname).(arc).vectorizedname{jj}=[basename '(' int2str(jj) ',:)'];
            end
        end

    case 'lagrangemultcc'
        cconstraintnum=retrieveimpulsemodelinformation(ocStruct,'inequalitycontrolconstraintnum');
        if ~cconstraintnum.value
            transStruct.(varname).name=[];
            transStruct.(varname).arrayname=[];
            transStruct.(varname).vectorizedname=[];
            return
        end
        basename=getbasicname(varname);
        transStruct.(varname).basicname=basename;
        for ii=1:cconstraintnum.value
            transStruct.(varname).name{ii}=[basename int2str(ii)];
            transStruct.(varname).arrayname{ii}=[basename '(' int2str(ii) ')'];
            transStruct.(varname).vectorizedname{ii}=[basename '(' int2str(ii) ',:)'];
        end
        
    case 'lagrangemultsc'
        sconstraintnum=retrieveimpulsemodelinformation(ocStruct,'inequalitystateconstraintnum');
        if ~sconstraintnum.value
            transStruct.(varname).name=[];
            transStruct.(varname).arrayname=[];
            transStruct.(varname).vectorizedname=[];
            return
        end
        basename=getbasicname(varname);
        transStruct.(varname).basicname=basename;
        for ii=1:sconstraintnum.value
            transStruct.(varname).name{ii}=[basename int2str(ii)];
            transStruct.(varname).arrayname{ii}=[basename '(' int2str(ii) ')'];
            transStruct.(varname).vectorizedname{ii}=[basename '(' int2str(ii) ',:)'];
        end
        
    case 'vectorizedsize'
        basename=getbasicname('dependent');
        transStruct.(varname).string=['size(' basename ',2)'];
        transStruct.(varname).arcdependent=0;

    case 'composedvar'
        transStruct.(varname).arcdependent=1;
        if ~isfield(transStruct,'state')
            transStruct=generatetransformationstruct4impulsemodel(ocStruct,transStruct,'state');
        end
        if ~isfield(transStruct,'statel')
            transStruct=generatetransformationstruct4impulsemodel(ocStruct,transStruct,'statel');
        end
        if ~isfield(transStruct,'stater')
            transStruct=generatetransformationstruct4impulsemodel(ocStruct,transStruct,'stater');
        end
        if ~isfield(transStruct,'costate')
            transStruct=generatetransformationstruct4impulsemodel(ocStruct,transStruct,'costate');
        end
        if ~isfield(transStruct,'costatel')
            transStruct=generatetransformationstruct4impulsemodel(ocStruct,transStruct,'costatel');
        end
        if ~isfield(transStruct,'costater')
            transStruct=generatetransformationstruct4impulsemodel(ocStruct,transStruct,'costater');
        end
        if ~isfield(transStruct,'control')
            transStruct=generatetransformationstruct4impulsemodel(ocStruct,transStruct,'control');
        end
        if ~isfield(transStruct,'impulsecontrol')
            transStruct=generatetransformationstruct4impulsemodel(ocStruct,transStruct,'impulsecontrol');
        end
        if ~isfield(transStruct,'lagrangemultcc')
            transStruct=generatetransformationstruct4impulsemodel(ocStruct,transStruct,'lagrangemultcc');
        end
        if ~isfield(transStruct,'lagrangemultsc')
            transStruct=generatetransformationstruct4impulsemodel(ocStruct,transStruct,'lagrangemultsc');
        end
        if ~isfield(transStruct,'exogenousfunction')
            transStruct=generatetransformationstruct4impulsemodel(ocStruct,transStruct,'exogenousfunction');
            transStruct=generatetransformationstruct4impulsemodel(ocStruct,transStruct,'exogenousfunctiondx');
        end
        for ii=1:ocStruct.arc.num
            if ii==1
                parametername=retrieveimpulsemodelinformation(ocStruct,'parametername');
                parameternum=retrieveimpulsemodelinformation(ocStruct,'parameternum');
                cconstraintnum=retrieveimpulsemodelinformation(ocStruct,'inequalitycontrolconstraintnum');
                sconstraintnum=retrieveimpulsemodelinformation(ocStruct,'inequalitystateconstraintnum');
                statename=retrieveimpulsemodelinformation(ocStruct,'statename');
                statenamel=retrieveimpulsemodelinformation(ocStruct,'statenamel');
                statenamer=retrieveimpulsemodelinformation(ocStruct,'statenamer');
                costatename=retrieveimpulsemodelinformation(ocStruct,'costatename');
                costatenamel=retrieveimpulsemodelinformation(ocStruct,'costatenamel');
                costatenamer=retrieveimpulsemodelinformation(ocStruct,'costatenamer');
                controlname=retrieveimpulsemodelinformation(ocStruct,'controlname');
                icontrolname=retrieveimpulsemodelinformation(ocStruct,'icontrolname');
                lagrangemultccname=retrieveimpulsemodelinformation(ocStruct,'lagrangemultipliercontrolname');
                exogenousfunctionnum=retrieveimpulsemodelinformation(ocStruct,'exogenousfunctionnum');
                exogenousfunctionname=retrieveimpulsemodelinformation(ocStruct,'exogenousfunctionnamewithargument');
                exogenousfunctionnamedx=retrieveimpulsemodelinformation(ocStruct,'exogenousfunctionnamedx');
                exogenousfunctionnamederivatived2x2=retrieveimpulsemodelinformation(ocStruct,'exogenousfunctionnamederivatived2x2');
                transStruct.(varname).parametername=parametername.value;
                for jj=1:parameternum.value
                    transStruct.(varname).vectorizedparametername{jj}=['repmat(' parametername.value{jj} ',1,size(' transStruct.state.basicname ',2))'];
                end
                transStruct.(varname).username=[statename.value statenamel.value statenamer.value costatename.value costatenamel.value costatenamer.value controlname.value icontrolname.value];
                transStruct.(varname).standardname=[transStruct.state.name transStruct.statel.name transStruct.stater.name transStruct.costate.name transStruct.costatel.name transStruct.costater.name transStruct.control.name transStruct.impulsecontrol.name];
                transStruct.(varname).standardarrayname=[transStruct.state.arrayname transStruct.statel.arrayname transStruct.stater.arrayname transStruct.costate.arrayname transStruct.costatel.arrayname transStruct.costater.arrayname transStruct.control.arrayname transStruct.impulsecontrol.arrayname];
                transStruct.(varname).standardvectorizedname=[transStruct.state.vectorizedname transStruct.statel.vectorizedname transStruct.stater.vectorizedname transStruct.costate.vectorizedname transStruct.costatel.vectorizedname transStruct.costater.vectorizedname transStruct.control.vectorizedname transStruct.impulsecontrol.vectorizedname];
                if cconstraintnum.value
                    transStruct.(varname).username=[transStruct.(varname).username lagrangemultccname.value];
                    transStruct.(varname).standardname=[transStruct.(varname).standardname transStruct.lagrangemultcc.name];
                    transStruct.(varname).standardarrayname=[transStruct.(varname).standardarrayname transStruct.lagrangemultcc.arrayname];
                    transStruct.(varname).standardvectorizedname=[transStruct.(varname).standardvectorizedname transStruct.lagrangemultcc.vectorizedname];
                end
                if sconstraintnum.value
                    transStruct.(varname).username=[transStruct.(varname).username lagrangemultscname.value];
                    transStruct.(varname).standardname=[transStruct.(varname).standardname transStruct.lagrangemultsc.name];
                    transStruct.(varname).standardarrayname=[transStruct.(varname).standardarrayname transStruct.lagrangemultsc.arrayname];
                    transStruct.(varname).standardvectorizedname=[transStruct.(varname).standardvectorizedname transStruct.lagrangemultsc.vectorizedname];
                end
                if exogenousfunctionnum.value
                    if ~isempty(exogenousfunctionnamedx.value)
                        transStruct.(varname).username=[transStruct.(varname).username exogenousfunctionname.value exogenousfunctionnamedx.value];
                        transStruct.(varname).standardname=[transStruct.(varname).standardname transStruct.exogenousfunction.name transStruct.exogenousfunctiondx.name];
                        transStruct.(varname).standardarrayname=[transStruct.(varname).standardarrayname transStruct.exogenousfunction.arrayname transStruct.exogenousfunctiondx.arrayname];
                        transStruct.(varname).standardvectorizedname=[transStruct.(varname).standardvectorizedname transStruct.exogenousfunction.vectorizedname transStruct.exogenousfunctiondx.vectorizedname];
                    else
                        transStruct.(varname).username=[transStruct.(varname).username exogenousfunctionname.value];
                        transStruct.(varname).standardname=[transStruct.(varname).standardname transStruct.exogenousfunction.name];
                        transStruct.(varname).standardarrayname=[transStruct.(varname).standardarrayname transStruct.exogenousfunction.arrayname];
                        transStruct.(varname).standardvectorizedname=[transStruct.(varname).standardvectorizedname transStruct.exogenousfunction.vectorizedname];
                    end
                end
            end
            arc=arcidentifier2field(ocStruct.arc.identifier{ii});
            if ~isempty(exogenousfunctionnamedx.value)
                transStruct.(varname).(arc).name=[transStruct.state.name transStruct.statel.name transStruct.stater.name transStruct.costate.name transStruct.costatel.name transStruct.costater.name transStruct.control.(arc).name transStruct.impulsecontrol.(arc).name transStruct.lagrangemultcc.name transStruct.lagrangemultsc.name transStruct.exogenousfunction.name transStruct.exogenousfunctiondx.name];
                transStruct.(varname).(arc).arrayname=[transStruct.state.arrayname transStruct.statel.arrayname transStruct.stater.arrayname transStruct.costate.arrayname transStruct.costatel.arrayname transStruct.costater.arrayname transStruct.control.(arc).arrayname transStruct.impulsecontrol.(arc).arrayname  transStruct.lagrangemultcc.arrayname  transStruct.lagrangemultsc.arrayname transStruct.exogenousfunction.arrayname transStruct.exogenousfunctiondx.arrayname];
                transStruct.(varname).(arc).vectorizedname=[transStruct.state.vectorizedname transStruct.statel.vectorizedname transStruct.stater.vectorizedname transStruct.costate.vectorizedname transStruct.costatel.vectorizedname transStruct.costater.vectorizedname transStruct.control.(arc).vectorizedname transStruct.impulsecontrol.(arc).vectorizedname  transStruct.lagrangemultcc.vectorizedname  transStruct.lagrangemultsc.vectorizedname transStruct.exogenousfunction.vectorizedname transStruct.exogenousfunctiondx.vectorizedname];
            else
                transStruct.(varname).(arc).name=[transStruct.state.name transStruct.statel.name transStruct.stater.name transStruct.costate.name transStruct.costatel.name transStruct.costater.name transStruct.control.(arc).name transStruct.lagrangemultcc.name transStruct.lagrangemultsc.name transStruct.exogenousfunction.name];
                transStruct.(varname).(arc).arrayname=[transStruct.state.arrayname transStruct.statel.arrayname transStruct.stater.arrayname transStruct.costate.arrayname transStruct.costatel.arrayname transStruct.costater.arrayname transStruct.control.(arc).arrayname  transStruct.lagrangemultcc.arrayname  transStruct.lagrangemultsc.arrayname transStruct.exogenousfunction.arrayname];
                transStruct.(varname).(arc).vectorizedname=[transStruct.state.vectorizedname transStruct.statel.vectorizedname transStruct.stater.vectorizedname transStruct.costate.vectorizedname transStruct.costatel.vectorizedname transStruct.costater.vectorizedname transStruct.control.(arc).vectorizedname  transStruct.lagrangemultcc.vectorizedname  transStruct.lagrangemultsc.vectorizedname transStruct.exogenousfunction.vectorizedname];
            end
        end
    case 'exogenousfunction'
        basename=getbasicname('exogenousfunction');
        exogenousfunctionnum=retrieveimpulsemodelinformation(ocStruct,'exogenousfunctionnum');
        transStruct.(varname).arcdependent=0;
        transStruct.(varname).basicname=basename;
        if exogenousfunctionnum.value
            for ii=1:exogenousfunctionnum.value
                transStruct.(varname).name{ii}=[basename int2str(ii)];
                transStruct.(varname).arrayname{ii}=[basename '(' int2str(ii) ')'];
                transStruct.(varname).vectorizedname{ii}=[basename '(' int2str(ii) ',:)'];
            end
        else
            transStruct.(varname).name=[];
            transStruct.(varname).arrayname=[];
            transStruct.(varname).vectorizedname=[];
        end

    case 'exogenousfunctiondx'
        basename=getbasicname('exogenousfunction');
        exogenousfunctionnum=retrieveimpulsemodelinformation(ocStruct,'exogenousfunctionnum');
        exogenousfunctionargument=retrieveimpulsemodelinformation(ocStruct,'exogenousfunctionargument');
        transStruct.(varname).arcdependent=0;
        transStruct.(varname).basicname=['D1' basename '_Dx'];
        if exogenousfunctionnum.value
            for ii=1:exogenousfunctionnum.value
                if ~isempty(exogenousfunctionargument.value{ii})
                    transStruct.(varname).name{ii}=['D1' basename '_Dx' int2str(ii)];
                    transStruct.(varname).arrayname{ii}=['D1' basename '_Dx(' int2str(ii) ')'];
                    transStruct.(varname).vectorizedname{ii}=['D1' basename '_Dx(' int2str(ii) ',:)'];
                else
                    transStruct.(varname).name{ii}='';
                    transStruct.(varname).arrayname{ii}='';
                    transStruct.(varname).vectorizedname{ii}='';
                end
            end
        else
            transStruct.(varname).name=[];
            transStruct.(varname).arrayname=[];
            transStruct.(varname).vectorizedname=[];
        end

end        