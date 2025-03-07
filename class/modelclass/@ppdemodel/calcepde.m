function epdeSol=calcepde(ppdeObj,coeff0,varargin)
%
% CALCEP calculates equilibria of an ocmodel.
%
% CALCEP(OCOBJ) calculates the equilibrium for the canonical system of the
% ocmodel OCOBJ using the symbolic toolbox.
%
% CALCEP(OCOBJ,X0) calculates the equilibrium numerically with starting
% values given at X0. If X0 is empty symbolic calculation is applied. For
% the numerical calculation 'FSOLVE' is used by default.
%
% CALCEP(OCOBJ,X0,ARCARG) equilibria are calculated for the canonical
% system corresponding to arc ARCARG. If ARCARG is empty or missing the
% equilibria for every arc are computed.
%
% CALCEP(OCOBJ,X0,ARCARG,OPT) in the ocoption structure OPT options for the
% numerical solver can be provided (category OPT.EQ). In the category
% GENERAL the equation sover can be changed
% opt.GENERAL.EquationSolver='fsolve' or ...
%
% EP = CALCEP(...) the output argument EP is a cell array of 'dynprimitive'
% instances.

global PPDEPRIMITIVE
opt=[]; % option for fsolve
epdeSol=[];
if isempty(ppdeObj)
    return
end
if nargin>=3
    opt=varargin{1};
end
if isempty(opt)
    opt=defaultocoptions;
end
funch=ellipticpde;
par=parametervalue(ppdeObj);
arcid=0;
equationsolver=getocoptions(opt,'GENERAL','EquationSolver');

switch equationsolver
    case 'fsolve'
        % test if optimization toolbox is registered
        if isempty(ver('optim'))
            ocmaterror('Optimization toolbox is not registered. Please provide another nonlinear eqaution solver.\n')
            return
        end
        [y,fval,exitflag,info,jacob]=fsolve(@zerofunction,coeff0,opt.EQ);
        ppdeTrrjStruct=struct(ppdetrajectory(PPDEPRIMITIVE.ppdePrimitive0));
        ppdeTrrjStruct.y=y;
        ppdeTrrjStruct.linearization=jacob;
        ppdeTrrjStruct.modelname=modelname(ppdeObj);
        ppdeTrrjStruct.modelparameter=par;
        epdeSol=ppdeprimitive(0,ppdetrajectory(ppdeTrrjStruct));
        %ppdeTrrjStruct.solver='fsolve';
end
    function [res J]=zerofunction(coeff)
        res=funch{1}(coeff,par,arcid);
        J=funch{3}(coeff,par,arcid);
    end
end
