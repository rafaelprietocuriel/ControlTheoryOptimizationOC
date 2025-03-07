function dataadaptation_bvp6c(tmesh)

global OCMATCONT OCBVP

nummesh=numel(tmesh);
arcposition=find(diff(tmesh)==0);
if OCMATCONT.HE.TIMEDDATA.nummesh==nummesh && (isempty(OCBVP.arcposition) && isempty(arcposition) || all(OCBVP.arcposition==arcposition))
    return
end
OCMATCONT.HE.TIMEDDATA.nummesh=nummesh;
% determine the number of arcs and the corresponding time mesh,
% switches of arcs are indicated by a doubling of the switch time
OCMATCONT.HE.TIMEDDATA.leftarcindex=[1 arcposition+1];
OCMATCONT.HE.TIMEDDATA.rightarcindex=[arcposition OCMATCONT.HE.TIMEDDATA.nummesh];
OCMATCONT.HE.TIMEDDATA.nummeshintv=OCMATCONT.HE.TIMEDDATA.rightarcindex-OCMATCONT.HE.TIMEDDATA.leftarcindex;
OCMATCONT.HE.DDATA.meshvalcoord=reshape(1:OCMATCONT.HE.TIMEDDATA.nummesh*OCBVP.neqn,OCBVP.neqn,OCMATCONT.HE.TIMEDDATA.nummesh);
nummeshcoord=OCMATCONT.HE.TIMEDDATA.nummesh*OCBVP.neqn;
OCMATCONT.HE.parametercoord=nummeshcoord+(1:OCMATCONT.HE.numparameter).';
OCMATCONT.HE.parametermcodcoord=nummeshcoord+(1:OCMATCONT.HE.numparameter-OCMATCONT.codimension).';
OCMATCONT.HE.contparametercoord=nummeshcoord+OCMATCONT.HE.numparameter;%-1+(1:OCMATCONT.codimension);
OCMATCONT.HE.numdvariables=OCMATCONT.HE.contparametercoord;
OCMATCONT.HE.numdvariablesmc=OCMATCONT.HE.numdvariables-OCMATCONT.codimension;
OCMATCONT.HE.coeffcoordmp=OCMATCONT.HE.numdvariables-OCMATCONT.HE.numparameter; % coordinates of coefficients without parameter values

OCBVP.arcposition=arcposition;
OCBVP.Lidx=OCMATCONT.HE.TIMEDDATA.leftarcindex;
OCBVP.Ridx=OCMATCONT.HE.TIMEDDATA.rightarcindex;
OCBVP.N=OCMATCONT.HE.TIMEDDATA.nummesh;
OCBVP.nN=OCBVP.numode*OCBVP.N;
OCBVP.Nint=OCMATCONT.HE.TIMEDDATA.nummeshintv;

