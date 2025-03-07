function ocObj=changespatialdiscretization(ocObj,N)
global ARCINFOVAR

if isempty(ocObj)
    return
end
basestatename=ocObj.Model.variable.state.name{1}(1:end-1);
basecostatename=ocObj.Model.variable.costate.name{1}(1:end-1);
basecontrolname=ocObj.Model.variable.control.name{1}(1:end-1);
for ii=0:N
    ocObj.Model.variable.state.name{ii+1}=[basestatename num2str(ii)];
    ocObj.Model.variable.costate.name{ii+1}=[basecostatename num2str(ii)];
    ocObj.Model.variable.control.name{ii+1}=[basecontrolname num2str(ii)];
end
ocObj.Model.variable.state.num=N+1;
ocObj.Model.variable.costate.num=N+1;
ocObj.Model.variable.control.num=N+1;
ocObj=changeparametervalue(ocObj,'N',N);

ARCINFOVAR.N=N;