function [tmesh,coeff,tangent,done,meshHistory,needGlobJac,EN]=meshadaptation(tmesh,coeff,tangent,PHI,done,meshHistory,needGlobJac,EN)
global OCMATCONT
OCMATCONT.prepare4meshadapt(tmesh,coeff);
% test if solution satisfies imposed conditions on the residual
res=OCMATCONT.residual(tmesh,coeff,PHI,OCMATCONT.ode);

maxres=max(res);
if OCMATCONT.verifyresidual(maxres)
    done = true;
elseif OCMATCONT.OPTIONS.meshadaptation % mesh adaptation is necessary
    % Detect mesh oscillations:  Was there a mesh with
    % the same number of nodes and a similar residual?
    residualReduction = abs(meshHistory((meshHistory(:,1) == OCMATCONT.HE.TIMEDDATA.nummesh),2) - maxres)/maxres;
    oscLikely = any( residualReduction < OCMATCONT.OPTIONS.residualreductionguard);

    % modify the mesh, interpolate the solution
    meshHistory(end+1,:) = [OCMATCONT.HE.TIMEDDATA.nummesh,maxres];
    [tmesh,coeff,tangent]=OCMATCONT.meshadaptation(tmesh,coeff,tangent,res,~oscLikely,OCMATCONT.ode);
    if OCMATCONT.HE.TIMEDDATA.nummesh>OCMATCONT.OPTIONS.maxgridnum
%         ocmatmsg('MATLAB:newtcorr4bvp:RelTolNotMet', ...
%             [ 'Unable to meet the tolerance without using more than %d '...
%             'mesh points. '],OCMATCONT.OPTIONS.maxgridnum);
        ocmatmsg([ 'Unable to meet the tolerance without using more than %d '...
            'mesh points. '],OCMATCONT.OPTIONS.maxgridnum);
        tmesh=[];
        coeff=[];
        tangent=[];
        return
    end
    needGlobJac = true;
    numdvariables=OCMATCONT.HE.numdvariables;
    EN=[];
    EN(numdvariables,1)=1;
    OCMATCONT.meshadapted=true;
else
%     ocmatmsg('MATLAB:newtcorr4bvp:RelTolNotMet', ...
%         'Residual tolerance %f is not met, but mesh is accepted with residual %f',OCMATCONT.OPTIONS.meshadaptreltol,maxres);
    ocmatmsg('Residual tolerance %f is not met, but mesh is accepted with residual %f\n',OCMATCONT.OPTIONS.meshadaptreltol,maxres);
    done = true;
end
