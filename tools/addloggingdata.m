function addloggingdata(loghdl,varargin)
global OCMATCONT

if ~ishandle(loghdl)
    ocmaterror('First Argument is not a handle')
    return
end
if nargin==1 || nargin==2 && isempty(varargin{1})
    LogBook.GeneralInfo.NewtonSolver.Method=OCMATCONT.newtonsolver;
    LogBook.GeneralInfo.NewtonSolver.AbsTol=OCMATCONT.OPTIONS.newtonabstol;
    LogBook.GeneralInfo.NewtonSolver.RelTol=OCMATCONT.OPTIONS.newtonreltol;
    LogBook.GeneralInfo.NewtonSolver.MaxNewtonIters=OCMATCONT.OPTIONS.maxnewtiter;
    LogBook.GeneralInfo.MeshAdaptation.MeshAdaptAbsTol=OCMATCONT.OPTIONS.meshadaptabstol;
    LogBook.GeneralInfo.MeshAdaptation.MeshAdaptRelTol=OCMATCONT.OPTIONS.meshadaptreltol;
    LogBook.GeneralInfo.MeshAdaptation.MeshAdaptMaxIter=OCMATCONT.OPTIONS.meshadaptmaxiter;
    LogBook.Continuation.NewtonSolver.FunctionEvaluation
    LogBook.Continuation.NewtonSolver.JacobianEvaluation
    LogBook.Continuation.NewtonSolver.NewtonIteration
    LogBook.Continuation.MeshAdaptation
    set(loghdl,'UserData',LogBook);
    return
end