function OCMATCONTC=commentedglobalocmatcont()
%
% returns structure where each field of the global variable for the ocmat
% continuation is explained


OCMATCONTC.DDATA.order='order of the ODEs of the canonical system';
OCMATCONT.DDATA.maxorder='maximal order';
OCMATCONT.DDATA.minorder='minimal order';
OCMATCONT.DDATA.sumorder='sum over all order';
OCMATCONTC.DDATA.numeq='number of the equations (ODEs,ADEQs)';
OCMATCONTC.DDATA.numeqcoord='1:OCMATCONTC.DDATA.numeq';
OCMATCONTC.DDATA.numcols='number of collocation points';
OCMATCONTC.DDATA.numcolscoord='1:OCMATCONTC.DDATA.numcols';
OCMATCONTC.DDATA.numparameter='number of free parameter values';
OCMATCONT.DDATA.multipointbvp='identifier if problem is forumalted as a multi-point BVP';
OCMATCONTC.factorial='values of the first 14 factorials';
OCMATCONT.DDATA.numdvariables='total number of variables for discretized problem';
% initialize global data for timegrid
OCMATCONTC.DDATA.mesh='mesh points (for the actual solution)';
OCMATCONTC.DDATA.nummesh='number of mesh points';
OCMATCONT.DDATA.t='time mesh';
OCMATCONTC.DDATA.dt='interval lengths of mesh [diff(x)]';
OCMATCONTC.DDATA.nummeshintv='number of mesh intervals';
OCMATCONTC.DDATA.cols='collocation points (for the actual solution)';
OCMATCONTC.DDATA.rho='normalized collocation points';
OCMATCONTC.DDATA.psi='Lagrange polynomials';
OCMATCONTC.DDATA.psival='Lagrange polynomials evaluated at collocation points';
OCMATCONTC.DDATA.intweights='Gaussian weights for integration';

OCMATCONTC.DDATA.gridvalcoord='coordinates for the values at the time mesh';
OCMATCONTC.DDATA.collvaldercoord='coordinates for the values of the derivatives at the collocation points';
OCMATCONTC.DDATA.parametercoord='coordinates for the parameter values';