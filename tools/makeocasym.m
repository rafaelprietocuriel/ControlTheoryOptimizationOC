function ocAsym=makeocasym(sol,varargin)

global OCMATCONT
problem='';
if nargin==2
    problem='extremal2ep';
end
if isempty(problem)
    if isempty(OCMATCONT)
        ocmaterror('Global variable of continuation proces has been cleared.')
    end
    problem=OCMATCONT.problem_func;
end
switch problem

    case 'extremal2ep'
        ocAsym=sol;
        funch=extremal2ep();
        pathname=funch{27}();
        load(fullfile(pathname,'SaveIntermediateResultsGlobalVariable4extremal2ep.mat'))
        ocAsym.limitset.y=MODELINFO.OCMATAE.saddlepoint;
        ocAsym.limitset.linearization=MODELINFO.OCMATAE.linearization;
        ocAsym=ocasymptotic(octrajectory(sol),dynprimitive(ocAsym.limitset));
end