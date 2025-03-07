function coeff0=initocmat_EPDE(ppdeObj,epdeSol0,varargin)
%
% INIT_EP_EP initialization for equilibrium continuation
%
% [X0,V0]= INIT_EP_EP(OCOBJ,DYNPRIM,AP)
clear global PPDEPRIMITIVE 
global PPDEPRIMITIVE 

if nargin==1
    epdeSol0=[];
end
if isempty(epdeSol0) % provide initialization function 
    epdeSol0=modelspecificfunc(ppdeObj,'InitPrimitiveSolution');
end
epdeSol0=ppdeprimitive(epdeSol0);

if isempty(ppdeObj) || isempty(epdeSol0)
    return
end
PPDEPRIMITIVE.modeltype=modeltype(ppdeObj);
PPDEPRIMITIVE.parametervalue=parametervalue(ppdeObj);
coeff0=dependentvar(epdeSol0);

switch PPDEPRIMITIVE.modeltype
    case 'ppdemodel'
        PPDEPRIMITIVE.pdefun=modelspecificfunc(ppdeObj,'CanonicalSystem');
        PPDEPRIMITIVE.pdejacobian=modelspecificfunc(ppdeObj,'CanonicalSystemJacobian');
        PPDEPRIMITIVE.parameterjacobian=modelspecificfunc(ppdeObj,'CanonicalSystemParameterJacobian');
        PPDEPRIMITIVE.points=points(epdeSol0);
        PPDEPRIMITIVE.edges=edges(epdeSol0);
        PPDEPRIMITIVE.triangles=triangles(epdeSol0);
        PPDEPRIMITIVE.boundarycondition=boundarycondition(epdeSol0);
        [M invM bcG K Kadv]=femoperator(epdeSol0);
        PPDEPRIMITIVE.femop.M=M;
        PPDEPRIMITIVE.femop.invM=full(invM);
        PPDEPRIMITIVE.femop.bcG=bcG;
        PPDEPRIMITIVE.femop.K=K;
        PPDEPRIMITIVE.femop.Kadv=Kadv;
        PPDEPRIMITIVE.numpoints=size(PPDEPRIMITIVE.points,2);
        idx=1:length(coeff0);

        N=length(idx)/PPDEPRIMITIVE.numpoints;
        PPDEPRIMITIVE.coeffidx=reshape(idx,PPDEPRIMITIVE.numpoints,N).';

end

modelfolder=getocmatfolder('userdata',modeltype(ppdeObj),modelname(ppdeObj));
PPDEPRIMITIVE.basicglobalvarfilename=fullocmatfile(modelfolder,'SaveIntermediateResults');
PPDEPRIMITIVE.basicresultfilename=fullocmatfile(modelfolder,'SaveIntermediateResultsGlobalVariable');
PPDEPRIMITIVE.parametervalue=parametervalue(ppdeObj);
PPDEPRIMITIVE.plotcoord=2;
PPDEPRIMITIVE.symjac=1;

PPDEPRIMITIVE.numdependentvar=2;
PPDEPRIMITIVE.ppdePrimitive0=epdeSol0;
