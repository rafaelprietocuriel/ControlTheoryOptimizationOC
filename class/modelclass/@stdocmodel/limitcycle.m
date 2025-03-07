function ocLC = limitcycle(ocObj)
%
% PERIODICSOLUTION returns the periodic solutions stored in OCOBJ
%
% OCEP=PERIODICSOLUTION(OCOBJ) returns a cell arry of periodic solutions
% stored in the  'Result' field 'LimitCycle' of OCOBJ.

% variable declaration
ocLC=dynprimitive();
if isempty(ocObj)
    return
end

ocResult=result(ocObj);
if isfield(ocResult,'LimitCycle')
    ocLC=ocResult.LimitCycle;
end
