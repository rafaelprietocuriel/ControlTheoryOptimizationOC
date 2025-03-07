function transStruct=generatetransformationstruct4standardmodelgrad(ocStruct,transStruct,varname,overwrite,varargin)
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
    case 'locstate'
        basename=getbasicname(varname);
        statenum=retrievemodelinformation(ocStruct,'statenum');
        transStruct.(varname).arcdependent=0;
        transStruct.(varname).basicname=basename;
        for ii=1:statenum.value
            transStruct.(varname).name{ii}=[basename int2str(ii)];
            transStruct.(varname).arrayname{ii}=[basename '(' int2str(ii) ')'];
            transStruct.(varname).vectorizedname{ii}=[basename '(' int2str(ii) ',:)'];
        end
    case 'loccostate'
        basename=getbasicname(varname);
        costatenum=retrievemodelinformation(ocStruct,'costatenum');
        transStruct.(varname).arcdependent=0;
        transStruct.(varname).basicname=basename;
        for ii=1:costatenum.value
            costatecoord=ii;
            transStruct.(varname).name{ii}=[basename int2str(costatecoord)];
            transStruct.(varname).arrayname{ii}=[basename '(' int2str(costatecoord) ')'];
            transStruct.(varname).vectorizedname{ii}=[basename '(' int2str(costatecoord) ',:)'];
        end

    case 'loccontrol'
        basename=getbasicname(varname);
        controlnum=retrievemodelinformation(ocStruct,'controlnum');
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

    case 'gradientloccontrol'
        basename=getbasicname(varname);
        controlnum=retrievemodelinformation(ocStruct,'controlnum');
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
    case 'vectorizedsize'
        basename=getbasicname('idependent');
        transStruct.(varname).string=['length(' basename ')'];
        transStruct.(varname).arcdependent=0;

    case 'locvar'
        transStruct.(varname).arcdependent=0;
        if ~isfield(transStruct,'locstate')
            transStruct=generatetransformationstruct4standardmodelgrad(ocStruct,transStruct,'locstate');
        end
        if ~isfield(transStruct,'loccostate')
            transStruct=generatetransformationstruct4standardmodelgrad(ocStruct,transStruct,'loccostate');
        end
        if ~isfield(transStruct,'loccontrol')
            transStruct=generatetransformationstruct4standardmodelgrad(ocStruct,transStruct,'loccontrol');
        end
        if ~isfield(transStruct,'gradientloccontrol')
            transStruct=generatetransformationstruct4standardmodelgrad(ocStruct,transStruct,'gradientloccontrol');
        end
        statename=retrievemodelinformation(ocStruct,'statename');
        costatename=retrievemodelinformation(ocStruct,'costatename');
        controlname=retrievemodelinformation(ocStruct,'controlname');
        if isempty(controlname.value)
            controlname.value='[]';
        end
        transStruct.(varname).username=[statename.value costatename.value controlname.value transStruct.gradientloccontrol.name];
        transStruct.(varname).standardname=[transStruct.locstate.name transStruct.loccostate.name transStruct.loccontrol.name transStruct.gradientloccontrol.name];
        transStruct.(varname).standardarrayname=[transStruct.locstate.arrayname transStruct.loccostate.arrayname transStruct.loccontrol.arrayname transStruct.gradientloccontrol.arrayname];
        transStruct.(varname).standardvectorizedname=[transStruct.locstate.vectorizedname transStruct.loccostate.vectorizedname transStruct.loccontrol.vectorizedname transStruct.gradientloccontrol.vectorizedname];
end
