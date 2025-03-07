function ocTrj=adapt2mdaesolver(ocTrj,tnew,collocationmeth,collocationnum,interpolationmeth)

opt=defaultocoptions;
if ~strcmp(solver(ocTrj),'dae')
    ocTrj=octrajectory2ocdae(ocTrj);
end

if isempty(ocTrj)
    return
end

if nargin<2
    collocationmeth=collocationmethod(ocTrj);
    if isempty(collocationmeth)
        collocationmeth=getocoptions(opt,'SBVPOC','CollocationMethod');
    end
end

if nargin<3
    collocationnum=collocationnumber(ocTrj);
    if isempty(collocationnum)
        collocationnum=getocoptions(opt,'SBVPOC','CollocationNumber');
    end
end

if nargin<4
    interpolationmeth=interpolationmethod(ocTrj);
    if isempty(interpolationmeth)
        interpolationmeth=getocoptions(opt,'SBVPOC','InterpolationMethod');
    end
end

sol.x=independentvar(ocTrj);
sol.y=dependentvar(ocTrj);
sol.data.coeff=rungekuttacoefficient(ocTrj);
sol.data.tfine=collocationmesh(ocTrj);
sol.parameters=freeparameter(ocTrj);

sol=adapt2mdaesolver(sol,tnew,daeorder(ocTrj),collocationmeth,collocationnum,interpolationmeth);

ocTrj.solverinfo.coeff=sol.data.coeff;
ocTrj.solverinfo.collocationmesh=sol.data.tfine;
ocTrj.solverinfo.collocationmethod=collocationmeth;
ocTrj.solverinfo.interpolationmethod=interpolationmeth;
ocTrj.solverinfo.collocationnumber=collocationnum;
