function dataadaptation_gbvp4c(tmesh)

global OCMATCONT OCBVP


nummesh=numel(tmesh);
arcposition=find(diff(tmesh)==0);
OCBVP.Lidxold=OCMATCONT.HE.TIMEDDATA.leftarcindex;
OCBVP.Ridxold=OCMATCONT.HE.TIMEDDATA.rightarcindex;
OCBVP.Nold=OCMATCONT.HE.TIMEDDATA.nummesh;
OCBVP.nNold=OCBVP.maxnumode*OCBVP.N;
OCBVP.Nintold=OCMATCONT.HE.TIMEDDATA.nummeshintv;
OCMATCONT.HE.DDATA.meshvalcoordold=OCMATCONT.HE.DDATA.meshvalcoord;
if OCMATCONT.HE.TIMEDDATA.nummesh==nummesh && (isempty(OCBVP.arcposition) && isempty(arcposition) || all(OCBVP.arcposition==arcposition))
    return
end

OCMATCONT.HE.TIMEDDATA.nummesh=nummesh;
% determine the number of arcs and the corresponding time mesh,
% switches of arcs are indicated by a doubling of the switch time
OCMATCONT.HE.TIMEDDATA.leftarcindex=[1 arcposition+1];
OCMATCONT.HE.TIMEDDATA.rightarcindex=[arcposition OCMATCONT.HE.TIMEDDATA.nummesh];
OCMATCONT.HE.TIMEDDATA.nummeshintv=OCMATCONT.HE.TIMEDDATA.rightarcindex-OCMATCONT.HE.TIMEDDATA.leftarcindex;


OCBVP.arcposition=arcposition;
OCBVP.Lidx=OCMATCONT.HE.TIMEDDATA.leftarcindex;
OCBVP.Ridx=OCMATCONT.HE.TIMEDDATA.rightarcindex;
OCBVP.N=OCMATCONT.HE.TIMEDDATA.nummesh;
OCBVP.nN=OCBVP.maxnumode*OCBVP.N;
OCBVP.Nint=zeros(1,size(OCBVP.Lidx,2));
for ii=1:size(OCBVP.Lidx,2)
    OCBVP.Nint(ii)=OCBVP.Ridx(ii)-OCBVP.Lidx(ii);
end
counter=0;
totidx=[];
for ii=1:size(OCBVP.Lidx,2)
    Y=zeros(OCBVP.maxnumode,OCBVP.Ridx(ii)-OCBVP.Lidx(ii)+1);
    Y(1:OCBVP.numode(ii),:)=1;
    idx=find(Y);
    totidx=[totidx;idx+counter];
    counter=counter+OCBVP.maxnumode*(OCBVP.Nint(ii)+1);
end
OCMATCONT.HE.DDATA.meshvalcoord=totidx;

%OCMATCONT.HE.DDATA.meshvalcoord=reshape(1:OCMATCONT.HE.TIMEDDATA.nummesh*OCBVP.maxnumode,OCBVP.maxnumode,OCMATCONT.HE.TIMEDDATA.nummesh);
nummeshcoord=length(totidx);
OCMATCONT.HE.parametercoord=nummeshcoord+(1:OCMATCONT.HE.numparameter).';
OCMATCONT.HE.parametermcodcoord=nummeshcoord+(1:OCMATCONT.HE.numparameter-OCMATCONT.codimension).';
OCMATCONT.HE.contparametercoord=nummeshcoord+OCMATCONT.HE.numparameter;%-1+(1:OCMATCONT.codimension);
OCMATCONT.HE.numdvariables=OCMATCONT.HE.contparametercoord;
OCMATCONT.HE.numdvariablesmc=OCMATCONT.HE.numdvariables-OCMATCONT.codimension;
OCMATCONT.HE.coeffcoordmp=OCMATCONT.HE.numdvariables-OCMATCONT.HE.numparameter; % coordinates of coefficients without parameter values
OCBVP.nN=sum(OCBVP.numode.*(OCBVP.Nint+1));
OCMATCONT.HE.ycoord=1:nummeshcoord;

