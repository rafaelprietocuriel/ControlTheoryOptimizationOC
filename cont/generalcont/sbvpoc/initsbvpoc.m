function solinit=initsbvpoc(ocObj,solObj,contclass,varargin)
% contclass: extremal, indifference, boundary
% extremal2ep, extremal2lc
solinit=[];

if ~issolutionclass(solObj)
    ocamterror('Second argument is not solutionclass.')
end
if isempty(ocObj) || ismepty(solObj)
    return
end

if strncmp(contclass,'extremal',8)
    initfunc=str2func(['initsbvpoc' contclass]);
end



solinit=initfunc(ocObj,solObj,contclass,varargin{:});