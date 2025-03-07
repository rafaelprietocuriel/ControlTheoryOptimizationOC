function ocPS = periodicsolution(ocObj)
%
% PERIODICSOLUTION returns the periodic solutions stored in OCOBJ
%
% OCEP=PERIODICSOLUTION(OCOBJ) returns a cell arry of periodic solutions
% stored in the  'Result' field 'PeriodicSolution' of OCOBJ.

% variable declaration
ocPS=dynprimitive();
if isempty(ocObj)
    return
end

ocResult=result(ocObj);
if isfield(ocResult,'PeriodicSolution')
    ocPS=ocResult.PeriodicSolution;
elseif isfield(ocResult,'LimitCycle')
    ocPS=ocResult.LimitCycle;
end
