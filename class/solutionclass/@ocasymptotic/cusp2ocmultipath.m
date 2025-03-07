function ocMP=cusp2ocmultipath(ocAsym,tau,v)

if nargin==2
    v=tangent(ocAsym);
end

if isempty(v)
    ocmatmessage('No tangent is provided.')
    ocMP=[];
    return
end

coeff=coefficient(ocAsym);
coeff1=coeff+tau*v;
coeff2=coeff-tau*v;
solveInfo=solverinfo(ocAsym);
ocAsym.octrajectory.solverinfo=[];
switch solveInfo.conttype
    case 'limitextremal'
        coeff1(end-1:end)=[];
        coeff2(end-1:end)=[];
end
ocMP{1}=ocAsym;
ocMP{1}.octrajectory.y=reshape(coeff1,size(ocAsym.octrajectory.y));
ocMP{2}=ocAsym;
ocMP{2}.octrajectory.y=reshape(coeff2,size(ocAsym.octrajectory.y));
ocMP=ocmultipath(ocMP);