function transStruct=generatetransformationstruct4standardmodel(ocStruct,transStruct,varname,overwrite,varargin)
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
        statenum=retrievemodelinformation(ocStruct,'statenum');
        transStruct.(varname).arcdependent=0;
        transStruct.(varname).basicname=basename;
        for ii=1:statenum.value
            transStruct.(varname).name{ii}=[basename int2str(ii)];
            transStruct.(varname).arrayname{ii}=[basename '(' int2str(ii) ')'];
            transStruct.(varname).vectorizedname{ii}=[basename '(' int2str(ii) ',:)'];
        end
    case 'costate'
        basename=getbasicname('dependent');
        costatenum=retrievemodelinformation(ocStruct,'costatenum');
        transStruct.(varname).arcdependent=0;
        transStruct.(varname).basicname=basename;
        for ii=1:costatenum.value
            costatecoord=ii+ocStruct.variable.state.num;
            transStruct.(varname).name{ii}=[basename int2str(costatecoord)];
            transStruct.(varname).arrayname{ii}=[basename '(' int2str(costatecoord) ')'];
            transStruct.(varname).vectorizedname{ii}=[basename '(' int2str(costatecoord) ',:)'];
        end
        if transform2statecontrolspace(ocStruct)
            arcnum=retrievemodelinformation(ocStruct,'arcnum');
            arcid=retrievemodelinformation(ocStruct,'arcidentifier');
            costatenum=retrievemodelinformation(ocStruct,'statenum');
            replacedbasename=getbasicname('costate');
            transStruct.(varname).arcdependent=1;
            for ii=1:arcnum.value
                arc=arcidentifier2field(arcid.value{ii});
                transform2statecontrolspaceflag=retrievemodelinformation(ocStruct,'transform2statecontrolspace',arcid.value{ii});
                if transform2statecontrolspaceflag.value
                    replacedcostateindex=retrievemodelinformation(ocStruct,'replacedcostateindex',arcid.value{ii});
                    for jj=1:costatenum.value
                        if any(jj==replacedcostateindex.value)
                            transStruct.(varname).(arc).name{jj}=[replacedbasename int2str(jj)];
                            transStruct.(varname).(arc).arrayname{jj}=[replacedbasename '(' int2str(jj) ')'];
                            transStruct.(varname).(arc).vectorizedname{jj}=[replacedbasename '(' int2str(jj) ',:)'];
                        else
                            costatecoord=jj+ocStruct.variable.state.num;
                            transStruct.(varname).(arc).name{jj}=[basename int2str(costatecoord)];
                            transStruct.(varname).(arc).arrayname{jj}=[basename '(' int2str(costatecoord) ')'];
                            transStruct.(varname).(arc).vectorizedname{jj}=[basename '(' int2str(costatecoord) ',:)'];
                        end
                    end
                else
                    transStruct.(varname).basicname=basename;
                    for jj=1:costatenum.value
                        costatecoord=ii+ocStruct.variable.state.num;
                        transStruct.(varname).(arc).name{jj}=[basename int2str(costatecoord)];
                        transStruct.(varname).(arc).arrayname{jj}=[basename '(' int2str(costatecoord) ')'];
                        transStruct.(varname).(arc).vectorizedname{jj}=[basename '(' int2str(costatecoord) ',:)'];
                    end
                end
            end
        end
    case 'variationstate'
        basename=getbasicname('dependent');
        variationnum=retrievemodelinformation(ocStruct,'variationparameternum');
        transStruct.(varname).arcdependent=0;
        transStruct.(varname).basicname=basename;
        if variationnum.value
            statenum=retrievemodelinformation(ocStruct,'statenum');
            costatenum=retrievemodelinformation(ocStruct,'costatenum');
            variationstatenum=retrievemodelinformation(ocStruct,'variationstatenum');
            shiftidx=statenum.value+costatenum.value;
            for ii=1:variationstatenum.value
                variationcoord=shiftidx+ii;
                transStruct.(varname).name{ii}=[basename int2str(variationcoord)];
                transStruct.(varname).arrayname{ii}=[basename '(' int2str(variationcoord) ')'];
                transStruct.(varname).vectorizedname{ii}=[basename '(' int2str(variationcoord) ',:)'];
            end
        else
            transStruct.(varname).name=[];
            transStruct.(varname).arrayname=[];
            transStruct.(varname).vectorizedname=[];
        end

    case 'variationcostate'
        basename=getbasicname('dependent');
        variationnum=retrievemodelinformation(ocStruct,'variationparameternum');
        transStruct.(varname).arcdependent=0;
        transStruct.(varname).basicname=basename;
        if variationnum.value
            statenum=retrievemodelinformation(ocStruct,'statenum');
            costatenum=retrievemodelinformation(ocStruct,'costatenum');
            variationstatenum=retrievemodelinformation(ocStruct,'variationstatenum');
            variationcostatenum=retrievemodelinformation(ocStruct,'variationcostatenum');
            transStruct.(varname).arcdependent=0;
            transStruct.(varname).basicname=basename;
            shiftidx=statenum.value+costatenum.value+variationstatenum.value;
            for ii=1:variationcostatenum.value
                variationcoord=shiftidx+ii;
                transStruct.(varname).name{ii}=[basename int2str(variationcoord)];
                transStruct.(varname).arrayname{ii}=[basename '(' int2str(variationcoord) ')'];
                transStruct.(varname).vectorizedname{ii}=[basename '(' int2str(variationcoord) ',:)'];
            end
        else
            transStruct.(varname).name=[];
            transStruct.(varname).arrayname=[];
            transStruct.(varname).vectorizedname=[];
        end

    case 'control'
        basename=getbasicname(varname);
        dynamicsbasename=getbasicname('dependent');
        controlnum=retrievemodelinformation(ocStruct,'controlnum');
        arcnum=retrievemodelinformation(ocStruct,'arcnum');
        statenum=retrievemodelinformation(ocStruct,'statenum');
        arcid=retrievemodelinformation(ocStruct,'arcidentifier');
        transStruct.(varname).arcdependent=1;
        transStruct.(varname).basicname=basename;
        if ~controlnum.value
            transStruct.(varname).name{1}='[]';
            transStruct.(varname).arrayname{1}='[]';
            transStruct.(varname).vectorizedname{1}='[]';
        end
        for ii=1:controlnum.value
            transStruct.(varname).name{ii}=[basename int2str(ii)];
            transStruct.(varname).arrayname{ii}=[basename '(' int2str(ii) ')'];
            transStruct.(varname).vectorizedname{ii}=[basename '(' int2str(ii) ',:)'];
        end
        totalimplicitvariableindex=retrievemodelinformation(ocStruct,'totalimplicitvariableindex');
        for ii=1:arcnum.value
            arc=arcidentifier2field(arcid.value{ii});
            implicitnonlinearcontrolindex=retrievemodelinformation(ocStruct,'implicitnonlinearcontrolindex',arcid.value{ii});
            implicitcounter=0;
            transform2statecontrolspaceflag=retrievemodelinformation(ocStruct,'transform2statecontrolspace',arcid.value{ii});
            if ~controlnum.value
                transStruct.(varname).(arc).name{1}='[]';
                transStruct.(varname).(arc).arrayname{1}='[]';
                transStruct.(varname).(arc).vectorizedname{1}='[]';
            end

            for jj=1:controlnum.value
                if any(implicitnonlinearcontrolindex.value==jj)
                    implicitcounter=implicitcounter+1;
                    controlcoord=2*statenum.value+implicitcounter;
                    transStruct.(varname).(arc).name{jj}=[dynamicsbasename int2str(controlcoord)];
                    transStruct.(varname).(arc).arrayname{jj}=[dynamicsbasename '(' int2str(controlcoord) ')'];
                    transStruct.(varname).(arc).vectorizedname{jj}=[dynamicsbasename '(' int2str(controlcoord) ',:)'];
                else
                    transStruct.(varname).(arc).name{jj}=[basename int2str(jj)];
                    transStruct.(varname).(arc).arrayname{jj}=[basename '(' int2str(jj) ')'];
                    transStruct.(varname).(arc).vectorizedname{jj}=[basename '(' int2str(jj) ',:)'];
                end
            end
            if transform2statecontrolspaceflag.value
                costatenum=retrievemodelinformation(ocStruct,'statenum');
                explicitcontroldynamicsindex=retrievemodelinformation(ocStruct,'explicitcontroldynamicsindex',arcid.value{ii});
                counter=0;
                for jj=1:controlnum.value
                    if any(explicitcontroldynamicsindex.value==jj)
                        counter=counter+1;
                        depcoord=costatenum.value+counter;
                        transStruct.(varname).(arc).name{jj}=[dynamicsbasename int2str(depcoord)];
                        transStruct.(varname).(arc).arrayname{jj}=[dynamicsbasename '(' int2str(depcoord) ')'];
                        transStruct.(varname).(arc).vectorizedname{jj}=[dynamicsbasename '(' int2str(depcoord) ',:)'];
                    else
                        transStruct.(varname).(arc).name{jj}=[basename int2str(jj)];
                        transStruct.(varname).(arc).arrayname{jj}=[basename '(' int2str(jj) ')'];
                        transStruct.(varname).(arc).vectorizedname{jj}=[basename '(' int2str(jj) ',:)'];
                    end
                end
            end
        end
        implicitcounter=0;
        for jj=1:controlnum.value
            if any(totalimplicitvariableindex.value==jj)
                implicitcounter=implicitcounter+1;
                controlcoord=2*statenum.value+implicitcounter;
                transStruct.(varname).stdarc.name{jj}=[dynamicsbasename int2str(controlcoord)];
                transStruct.(varname).stdarc.arrayname{jj}=[dynamicsbasename '(' int2str(controlcoord) ')'];
                transStruct.(varname).stdarc.vectorizedname{jj}=[dynamicsbasename '(' int2str(controlcoord) ',:)'];
            else
                transStruct.(varname).stdarc.name{jj}=[basename int2str(jj)];
                transStruct.(varname).stdarc.arrayname{jj}=[basename '(' int2str(jj) ')'];
                transStruct.(varname).stdarc.vectorizedname{jj}=[basename '(' int2str(jj) ',:)'];
            end
        end

    case 'lagrangemultcc'
        cconstraintnum=retrievemodelinformation(ocStruct,'inequalitycontrolconstraintnum');
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

    case 'controldynamics'
        basename=getbasicname('control');
        basename=['D' basename 'Dt'];
        controlnum=retrievemodelinformation(ocStruct,'controlnum');
        transStruct.(varname).arcdependent=0;
        transStruct.(varname).basicname=basename;
        for ii=1:controlnum.value
            transStruct.(varname).name{ii}=[basename int2str(ii)];
            transStruct.(varname).arrayname{ii}=[basename '(' int2str(ii) ')'];
            transStruct.(varname).vectorizedname{ii}=[basename '(' int2str(ii) ',:)'];
        end
        cconstraintnum=retrievemodelinformation(ocStruct,'inequalitycontrolconstraintnum');
        if cconstraintnum.value
            implicitnonlinearcontrolindex=retrievemodelinformation(ocStruct,'totalimplicitvariableindex');
            if implicitnonlinearcontrolindex.value
                controlnum=retrievemodelinformation(ocStruct,'controlnum');
                implicitcounter=controlnum.value;
                for jj=1:cconstraintnum.value
                    implicitcounter=implicitcounter+1;
                    transStruct.(varname).name{controlnum.value+jj}=[basename int2str(implicitcounter)];
                    transStruct.(varname).arrayname{controlnum.value+jj}=[basename '(' int2str(implicitcounter) ')'];
                    transStruct.(varname).vectorizedname{controlnum.value+jj}=[basename '(' int2str(implicitcounter) ',:)'];
                end
            end
        end
    case 'lagrangemultsc'
        sconstraintnum=retrievemodelinformation(ocStruct,'inequalitystateconstraintnum');
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
        basename=getbasicname('independent');
        transStruct.(varname).string=['length(' basename ')'];
        transStruct.(varname).arcdependent=0;

    case 'exogenousstate'
        basename=getbasicname('exogenousstate');
        exogenousdynamicsnum=retrievemodelinformation(ocStruct,'exogenousdynamicsnum');
        statenum=retrievemodelinformation(ocStruct,'statenum');
        costatenum=retrievemodelinformation(ocStruct,'costatenum');
        transStruct.(varname).arcdependent=0;
        transStruct.(varname).basicname=basename;
        if exogenousdynamicsnum.value
            for ii=1:exogenousdynamicsnum.value
                coord=statenum.value+costatenum.value+ii;
                transStruct.(varname).name{ii}=[basename int2str(coord)];
                transStruct.(varname).arrayname{ii}=[basename '(' int2str(coord) ')'];
                transStruct.(varname).vectorizedname{ii}=[basename '(' int2str(coord) ',:)'];
            end
        else
            transStruct.(varname).name=[];
            transStruct.(varname).arrayname=[];
            transStruct.(varname).vectorizedname=[];
        end

    case 'exogenousfunction'
        basename=getbasicname('exogenousfunction');
        exogenousfunctionnum=retrievemodelinformation(ocStruct,'exogenousfunctionnum');
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
        exogenousfunctionnum=retrievemodelinformation(ocStruct,'exogenousfunctionnum');
        exogenousfunctionargument=retrievemodelinformation(ocStruct,'exogenousfunctionargument');
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

    case 'nonsmoothfunction'
        basename=getbasicname('nonsmoothfunction');
        nonsmoothfunctionnum=retrievemodelinformation(ocStruct,'nonsmoothfunctionnum');
        transStruct.(varname).arcdependent=0;
        transStruct.(varname).basicname=basename;
        if nonsmoothfunctionnum.value
            for ii=1:nonsmoothfunctionnum.value
                transStruct.(varname).name{ii}=[basename int2str(ii)];
                transStruct.(varname).arrayname{ii}=[basename '(' int2str(ii) ')'];
                transStruct.(varname).vectorizedname{ii}=[basename '(' int2str(ii) ',:)'];
            end
        else
            transStruct.(varname).name=[];
            transStruct.(varname).arrayname=[];
            transStruct.(varname).vectorizedname=[];
        end

    case 'nonsmoothfunctiondx'
        basename=getbasicname('nonsmoothfunction');
        nonsmoothfunctionnum=retrievemodelinformation(ocStruct,'nonsmoothfunctionnum');
        nonsmoothfunctionargument=retrievemodelinformation(ocStruct,'allnonsmoothfunctionargument');
        statenum=retrievemodelinformation(ocStruct,'statenum');
        transStruct.(varname).arcdependent=0;
        transStruct.(varname).basicname=['D1' basename '_Dx'];
        if nonsmoothfunctionnum.value
            for ii=1:nonsmoothfunctionnum.value
                if ~isempty(nonsmoothfunctionargument.value{ii})
                    for jj=1:statenum.value
                        transStruct.(varname).name{ii}{jj}=['D1' basename '_Dx' int2str(jj)];
                        transStruct.(varname).arrayname{ii}{jj}=['D1' basename '_Dx(' int2str(jj) ')'];
                        transStruct.(varname).vectorizedname{ii}{jj}=['D1' basename '_Dx(' int2str(jj) ',:,' int2str(ii) ')'];
                    end
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
    case 'composedvar'
        transStruct.(varname).arcdependent=1;
        if ~isfield(transStruct,'independent')
            transStruct=generatetransformationstruct4standardmodel(ocStruct,transStruct,'independent');
        end
        if ~isfield(transStruct,'state')
            transStruct=generatetransformationstruct4standardmodel(ocStruct,transStruct,'state');
        end
        if ~isfield(transStruct,'costate')
            transStruct=generatetransformationstruct4standardmodel(ocStruct,transStruct,'costate');
        end
        if ~isfield(transStruct,'control')
            transStruct=generatetransformationstruct4standardmodel(ocStruct,transStruct,'control');
        end
        if ~isfield(transStruct,'lagrangemultcc')
            transStruct=generatetransformationstruct4standardmodel(ocStruct,transStruct,'lagrangemultcc');
        end
        if ~isfield(transStruct,'lagrangemultsc')
            transStruct=generatetransformationstruct4standardmodel(ocStruct,transStruct,'lagrangemultsc');
        end
        if ~isfield(transStruct,'exogenousstate')
            transStruct=generatetransformationstruct4standardmodel(ocStruct,transStruct,'exogenousstate');
        end
        if ~isfield(transStruct,'exogenousfunction')
            transStruct=generatetransformationstruct4standardmodel(ocStruct,transStruct,'exogenousfunction');
            transStruct=generatetransformationstruct4standardmodel(ocStruct,transStruct,'exogenousfunctiondx');
        end
        if ~isfield(transStruct,'nonsmoothfunction')
            transStruct=generatetransformationstruct4standardmodel(ocStruct,transStruct,'nonsmoothfunction');
            transStruct=generatetransformationstruct4standardmodel(ocStruct,transStruct,'nonsmoothfunctiondx');
        end
        if ~isfield(transStruct,'variationstate')
            transStruct=generatetransformationstruct4standardmodel(ocStruct,transStruct,'variationstate');
        end
        if ~isfield(transStruct,'variationcostate')
            transStruct=generatetransformationstruct4standardmodel(ocStruct,transStruct,'variationcostate');
        end
        if implicitcontrols(ocStruct)
            if ~isfield(transStruct,'controldynamics')
                transStruct=generatetransformationstruct4standardmodel(ocStruct,transStruct,'controldynamics');
            end
        end

        for ii=1:ocStruct.arc.num
            if ii==1
                parametername=retrievemodelinformation(ocStruct,'parametername');
                parameternum=retrievemodelinformation(ocStruct,'parameternum');
                cconstraintnum=retrievemodelinformation(ocStruct,'inequalitycontrolconstraintnum');
                sconstraintnum=retrievemodelinformation(ocStruct,'inequalitystateconstraintnum');
                statename=retrievemodelinformation(ocStruct,'statename');
                costatename=retrievemodelinformation(ocStruct,'costatename');
                controlname=retrievemodelinformation(ocStruct,'controlname');
                lagrangemultccname=retrievemodelinformation(ocStruct,'lagrangemultipliercontrolname');
                lagrangemultscname=retrievemodelinformation(ocStruct,'lagrangemultiplierstatename');
                nonsmoothfunctionnum=retrievemodelinformation(ocStruct,'nonsmoothfunctionnum');
                nonsmoothfunctionname=retrievemodelinformation(ocStruct,'allnonsmoothfunctionnamewithbraces');
                nonsmoothfunctionnamedx=retrievemodelinformation(ocStruct,'allnonsmoothfunctionnamedx');
                exogenousdynamicsnum=retrievemodelinformation(ocStruct,'exogenousdynamicsnum');
                exogenousstatename=retrievemodelinformation(ocStruct,'exogenousstatename');
                exogenousfunctionnum=retrievemodelinformation(ocStruct,'exogenousfunctionnum');
                exogenousfunctionname=retrievemodelinformation(ocStruct,'exogenousfunctionnamewithargument');
                exogenousfunctionnamedx=retrievemodelinformation(ocStruct,'exogenousfunctionnamedx');
                generalizedcontroldynamicsvariable=retrievemodelinformation(ocStruct,'generalizedcontroldynamicsvariable');
                implicitnonlinearcontrolindex=retrievemodelinformation(ocStruct,'totalimplicitvariableindex');
                variationnum=retrievemodelinformation(ocStruct,'variationparameternum');
                if variationnum.value
                    variationstatename=retrievemodelinformation(ocStruct,'variationstatename');
                    variationcostatename=retrievemodelinformation(ocStruct,'variationcostatename');
                end
                transStruct.(varname).parametername=parametername.value;
                for jj=1:parameternum.value
                    transStruct.(varname).vectorizedparametername{jj}=['repmat(' parametername.value{jj} ',1,length(' transStruct.independent.name '))'];
                end
                if isempty(controlname.value)
                    controlname.value='[]';
                end
                transStruct.(varname).username=[statename.value costatename.value controlname.value];
                transStruct.(varname).standardname=[transStruct.state.name transStruct.costate.name transStruct.control.name];
                transStruct.(varname).standardarrayname=[transStruct.state.arrayname transStruct.costate.arrayname transStruct.control.arrayname];
                transStruct.(varname).standardvectorizedname=[transStruct.state.vectorizedname transStruct.costate.vectorizedname transStruct.control.vectorizedname];
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
                if exogenousdynamicsnum.value
                    transStruct.(varname).username=[transStruct.(varname).username exogenousstatename.value];
                    transStruct.(varname).standardname=[transStruct.(varname).standardname transStruct.exogenousstate.name];
                    transStruct.(varname).standardarrayname=[transStruct.(varname).standardarrayname transStruct.exogenousstate.arrayname];
                    transStruct.(varname).standardvectorizedname=[transStruct.(varname).standardvectorizedname transStruct.exogenousstate.vectorizedname];
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
                if variationnum.value
                    transStruct.(varname).username=[transStruct.(varname).username variationstatename.value variationcostatename.value];
                    transStruct.(varname).standardname=[transStruct.(varname).standardname transStruct.variationstate.name transStruct.variationcostate.name];
                    transStruct.(varname).standardarrayname=[transStruct.(varname).standardarrayname transStruct.variationstate.arrayname transStruct.variationcostate.arrayname];
                    transStruct.(varname).standardvectorizedname=[transStruct.(varname).standardvectorizedname transStruct.variationstate.vectorizedname transStruct.variationcostate.vectorizedname];
                end
                if nonsmoothfunctionnum.value
                    if ~isempty(nonsmoothfunctionnamedx.value)
                        transStruct.(varname).username=[transStruct.(varname).username nonsmoothfunctionname.value [nonsmoothfunctionnamedx.value{:}]];
                        transStruct.(varname).standardname=[transStruct.(varname).standardname transStruct.nonsmoothfunction.name [transStruct.nonsmoothfunctiondx.name{:}]];
                        transStruct.(varname).standardarrayname=[transStruct.(varname).standardarrayname transStruct.nonsmoothfunction.arrayname [transStruct.nonsmoothfunctiondx.arrayname{:}]];
                        transStruct.(varname).standardvectorizedname=[transStruct.(varname).standardvectorizedname transStruct.nonsmoothfunction.vectorizedname [transStruct.nonsmoothfunctiondx.vectorizedname{:}]];
                    else
                        transStruct.(varname).username=[transStruct.(varname).username nonsmoothfunctionname.value];
                        transStruct.(varname).standardname=[transStruct.(varname).standardname transStruct.nonsmoothfunction.name];
                        transStruct.(varname).standardarrayname=[transStruct.(varname).standardarrayname transStruct.nonsmoothfunction.arrayname];
                        transStruct.(varname).standardvectorizedname=[transStruct.(varname).standardvectorizedname transStruct.nonsmoothfunction.vectorizedname];
                    end
                end
                if ~isempty(implicitnonlinearcontrolindex.value)
                    transStruct.(varname).username=[transStruct.(varname).username generalizedcontroldynamicsvariable.value];
                    transStruct.(varname).standardname=[transStruct.(varname).standardname transStruct.controldynamics.name];
                    transStruct.(varname).standardarrayname=[transStruct.(varname).standardarrayname transStruct.controldynamics.arrayname];
                    transStruct.(varname).standardvectorizedname=[transStruct.(varname).standardvectorizedname transStruct.controldynamics.vectorizedname];
                end
            end
            arc=arcidentifier2field(ocStruct.arc.identifier{ii});
            if ~isempty(exogenousfunctionnamedx.value)
                transStruct.(varname).(arc).name=[transStruct.state.name transStruct.costate.name transStruct.control.(arc).name transStruct.lagrangemultcc.name transStruct.lagrangemultsc.name transStruct.exogenousfunction.name transStruct.exogenousfunctiondx.name transStruct.variationstate.name transStruct.variationcostate.name];
                transStruct.(varname).(arc).arrayname=[transStruct.state.arrayname transStruct.costate.arrayname transStruct.control.(arc).arrayname  transStruct.lagrangemultcc.arrayname  transStruct.lagrangemultsc.arrayname transStruct.exogenousfunction.arrayname transStruct.exogenousfunctiondx.arrayname transStruct.variationstate.arrayname transStruct.variationcostate.arrayname];
                transStruct.(varname).(arc).vectorizedname=[transStruct.state.vectorizedname transStruct.costate.vectorizednamee transStruct.control.(arc).vectorizedname  transStruct.lagrangemultcc.vectorizedname  transStruct.lagrangemultsc.vectorizedname transStruct.exogenousfunction.vectorizedname transStruct.exogenousfunctiondx.vectorizedname transStruct.variationstate.vectorizedname transStruct.variationcostate.vectorizedname];
            else
                transStruct.(varname).(arc).name=[transStruct.state.name transStruct.costate.name transStruct.control.(arc).name transStruct.lagrangemultcc.name transStruct.lagrangemultsc.name transStruct.exogenousfunction.name transStruct.variationstate.name transStruct.variationcostate.name transStruct.nonsmoothfunction.name];
                transStruct.(varname).(arc).arrayname=[transStruct.state.arrayname transStruct.costate.arrayname transStruct.control.(arc).arrayname  transStruct.lagrangemultcc.arrayname  transStruct.lagrangemultsc.arrayname transStruct.exogenousfunction.arrayname transStruct.variationstate.arrayname transStruct.variationcostate.arrayname transStruct.nonsmoothfunction.arrayname];
                transStruct.(varname).(arc).vectorizedname=[transStruct.state.vectorizedname transStruct.costate.vectorizedname transStruct.control.(arc).vectorizedname  transStruct.lagrangemultcc.vectorizedname  transStruct.lagrangemultsc.vectorizedname transStruct.exogenousfunction.vectorizedname transStruct.variationstate.vectorizedname transStruct.variationcostate.vectorizedname transStruct.nonsmoothfunction.vectorizedname];
            end
            if exogenousdynamicsnum.value
                transStruct.(varname).(arc).name=[transStruct.(varname).(arc).name transStruct.exogenousstate.name];
                transStruct.(varname).(arc).arrayname=[transStruct.(varname).(arc).arrayname transStruct.exogenousstate.arrayname];
                transStruct.(varname).(arc).vectorizedname=[transStruct.(varname).(arc).vectorizedname transStruct.exogenousstate.vectorizedname];
            end
            if ~isempty(implicitnonlinearcontrolindex.value)
                transStruct.(varname).(arc).name=[transStruct.(varname).(arc).name transStruct.controldynamics.name];
                transStruct.(varname).(arc).arrayname=[transStruct.(varname).(arc).arrayname transStruct.controldynamics.arrayname];
                transStruct.(varname).(arc).vectorizedname=[transStruct.(varname).(arc).vectorizedname transStruct.controldynamics.vectorizedname];
            end

        end
        if ~isempty(implicitnonlinearcontrolindex.value)
            transStruct.(varname).stdarc.name=[transStruct.state.name transStruct.costate.name transStruct.control.stdarc.name transStruct.lagrangemultcc.name transStruct.lagrangemultsc.name transStruct.exogenousfunction.name transStruct.exogenousfunctiondx.name];
            transStruct.(varname).stdarc.arrayname=[transStruct.state.arrayname transStruct.costate.arrayname transStruct.control.stdarc.arrayname  transStruct.lagrangemultcc.arrayname  transStruct.lagrangemultsc.arrayname transStruct.exogenousfunction.arrayname transStruct.exogenousfunctiondx.arrayname];
            transStruct.(varname).stdarc.vectorizedname=[transStruct.state.vectorizedname transStruct.costate.vectorizedname transStruct.control.stdarc.vectorizedname  transStruct.lagrangemultcc.vectorizedname  transStruct.lagrangemultsc.vectorizedname transStruct.exogenousfunction.vectorizedname transStruct.exogenousfunctiondx.vectorizedname];
        end
end