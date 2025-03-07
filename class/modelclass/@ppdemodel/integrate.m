function out=integrate(ppdeObj,functionname,solObj,integrationvariable,varargin)
%
out=[];
if isempty(ppdeObj) || isempty(solObj)
    return
end

spacearg=spaceargument(ppdeObj);
timearg=timeargument(ppdeObj);
switch integrationvariable
    case spacearg
        spacedim=spacedimension(ppdeObj);
        switch spacedim
            case 1
                femdat=femdata(solObj);
                z=femdat.grid.x;
                dz=diff(z);
                dz=dz(:);
                zmid=(z(1:end-1)+z(2:end))*0.5;
                if ischar(functionname)
                    SolObjMid=spacedeval(solObj,zmid,ppdeObj);
                    fval=eval([functionname '(ppdeObj,SolObjMid,1);']);
                elseif isnumeric(functionname)
                    fval=dependentvar(solObj);
                    fval=interp1(z,fval(functionname,:),zmid,'linear');
                else
                end
            otherwise
                return
        end
    case timearg
end
out=sum(fval.*repmat(dz,1,size(fval,2)));