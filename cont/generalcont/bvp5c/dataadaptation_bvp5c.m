function dataadaptation_bvp5c(tmesh)

global OCMATCONT OCBVP

OCMATCONT.HE.TIMEDDATA.nummesh=numel(tmesh);
% determine the number of arcs and the corresponding time mesh,
% switches of arcs are indicated by a doubling of the switch time
OCMATCONT.HE.TIMEDDATA.leftarcindex=[1 1];
OCMATCONT.HE.TIMEDDATA.rightarcindex=[1 OCMATCONT.HE.TIMEDDATA.nummesh];
OCMATCONT.HE.TIMEDDATA.nummeshintv=OCMATCONT.HE.TIMEDDATA.rightarcindex-OCMATCONT.HE.TIMEDDATA.leftarcindex;
%OCMATCONT.HE.numarc=1;
OCMATCONT.HE.DDATA.meshvalcoord=reshape(1:OCMATCONT.HE.TIMEDDATA.nummesh*OCBVP.neqn,OCBVP.neqn,OCMATCONT.HE.TIMEDDATA.nummesh);
OCMATCONT.HE.parametercoord=OCBVP.n+(1:OCMATCONT.HE.numparameter).';
OCMATCONT.HE.contparametercoord=OCBVP.neqn-1+(1:OCMATCONT.codimension);
OCMATCONT.HE.numdvariables=numel(OCMATCONT.HE.DDATA.meshvalcoord);
OCMATCONT.HE.numdvariablesmc=OCMATCONT.HE.numdvariables-OCMATCONT.codimension;
OCMATCONT.HE.coeffcoordmp=OCMATCONT.HE.numdvariables-OCMATCONT.HE.numparameter; % coordinates of coefficients without parameter values

OCBVP.N=OCMATCONT.HE.TIMEDDATA.nummesh;
OCBVP.nN=OCBVP.n*OCBVP.N;
