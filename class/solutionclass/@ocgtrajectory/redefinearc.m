function ocgTrj=redefinearc(ocgTrj,newposition,arcid,varargin)

ocObj=stdocmodel(modelname(ocgTrj),[],[],0);
par=modelparameter(ocgTrj);
ocObj=changeparametervalue(ocObj,par);
staten=statenum(ocObj);
ctrln=controlnum(ocObj);
maxdim=2*staten+ctrln;
ocTrj=struct(ocgTrj.octrajectory);
arcpos=arcposition(ocgTrj);
arcn=arcnum(ocgTrj);

y=zeros(maxdim,arcpos(2,end));
ctrl=control(ocgTrj,ocObj);
for ii=1:arcn
    y(1:maxdim,arcpos(1,ii):arcpos(2,ii))=[ocTrj.y{ii}(1:2*staten,:);ctrl{ii}];
end
ocTrj.y=y;
ocTrj=octrajectory(redefinearc(ocTrj,newposition,arcid,varargin{:}));
arcn=arcnum(ocTrj);
arcpos=arcposition(ocTrj);
arcarg=arcargument(ocTrj);
y=ocTrj.y;
ocTrj.y=[];
odenum=zeros(1,arcn);
for ii=1:arcn
    icoord=implicitcontrolcoordinate(ocObj,arcarg(ii));
    if ~isempty(icoord)
        ocTrj.y{ii}=y([1:2*staten 2*staten+icoord],arcpos(1,ii):arcpos(2,ii));
    else
        ocTrj.y{ii}=y(1:2*staten,arcpos(1,ii):arcpos(2,ii));
    end
    odenum(ii)=canonicalsystemdimension(ocObj,arcarg(ii));
end
ocgTrj=struct(ocTrj);
ocgTrj.odenum=odenum;
ocgTrj=ocgtrajectory(ocgTrj);