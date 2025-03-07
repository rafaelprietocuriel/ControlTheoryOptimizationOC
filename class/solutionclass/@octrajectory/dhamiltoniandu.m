function varargout=dhamiltoniandu(ocTrj,connectflag)
%
%
if nargin==1
    connectflag=0;
end

ocObj=loadmodel(ocTrj);
if isempty(ocObj)
    varargout{1}=[];
    return
end
par=modelparameter(ocTrj);
ocObj=changeparametervalue(ocObj,par);

arcarg=arcargument(ocTrj);
indepvar=time(ocObj,ocTrj,1);
depvar=dependentvar(ocTrj);
arcpos=arcposition(ocTrj);
arcn=arcnum(ocTrj);
ctrlnum=controlnum(ocObj);
cstrnum=inequalitycontrolconstraintnum(ocObj);


for ii=1:arcn
    actcstrnum=sum(constraintcombinationindex(ocObj,arcarg));
    arcp=arcpos(1,ii):arcpos(2,ii);
    dHdu=feval(ocObj,'DHamiltonianDU',indepvar(arcpos(1,ii):arcpos(2,ii)),depvar(:,arcpos(1,ii):arcpos(2,ii)),par,arcarg(ii));
    if connectflag
        varargout{1}(1:(ctrlnum+cstrnum),arcp)=zeros(ctrlnum+cstrnum,arcpos(2,ii)-arcpos(1,ii)+1);
        varargout{1}(1:(ctrlnum+actcstrnum),arcp)=dHdu;
    else
        varargout{ii}=dHdu;
    end
end
