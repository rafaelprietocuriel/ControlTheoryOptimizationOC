function [depvarnew,ctrl_locnew,arcidnew]=ocshootloc(t,depvar,ctrl_loc,par,arcid)
% solves a (local) optimal control problem using the shooting method
%
% for the shooting method the control is handled like an exogenously given
% time function

global OCSHOOTCONT OCGRADSOL

% solve canonical system with exogenously given control ctrl_loc
depvarnew=OCSHOOTCONT.shootingodesolver(OCGRADSOL.canonicalsystemdynamics,t,depvar(:,1),ctrl_loc,par);

if OCSHOOTCONT.cconstraint
    ctr=0;
    while 1
        ctr=ctr+1;
        [ctrl_locnew,arcidnew]=correctcontrol(t,depvarnew,par);
        depvarnew=OCSHOOTCONT.shootingodesolver(OCGRADSOL.canonicalsystemdynamics,t,depvarnew(:,1),ctrl_locnew,par);
        admval=isadmissible(t,depvarnew,ctrl_locnew,par,arcidnew);
%         if ctr>2
%             ctr
%         end
        if admval>-1e-4 || ctr>10
            if ctr==10
                [ctrl_locnew,arcidnew]=correctcontrol(t,depvarnew,par);
            end
            break
        end
    end
else
    [ctrl_locnew,arcidnew]=correctcontrol(t,depvarnew,par);
end

function [ctrl_loc,arcid]=correctcontrol(t,depvar,par)
global OCSHOOTCONT OCGRADSOL

switch OCSHOOTCONT.OPTIONS.correctcontroltype
    case 'prior' % the Hamiltonian maximizing condition can be solved ahead
        
        % return extended Hamiltonian value for each arcid
        L=OCGRADSOL.extendedhamiltonian(t,depvar,par); % L~R^NxR^nt, N number of constraint combinations, nt time-grid number
        L(abs(imag(L))>0)=NaN;
        
        % return Lagrangian multipliers for each arcid
        lm=OCGRADSOL.lagrangemultiplier(t,depvar,par); %lm~R^NxR^ntxR^an % an number of arcids
        cstr=OCGRADSOL.constraint(t,depvar,par); %lm~R^NxR^ntxR^an % an number of arcids

        % exclude values that are either imaginary or are negative
        tmp=squeeze((sum(imag(abs(lm))>0,1)+sum(lm<0,1)+sum(isnan(lm),1))).'+squeeze((sum(imag(abs(cstr))>0,1)+sum(cstr<0,1)+sum(isnan(cstr),1))).';
        

        L(tmp>0)=NaN;

        % return maximized extended hamiltonian and the according arcid
        [Lmax,arcid]=max(L,[],1);
        arcid=arcid-1;
        if numel(~isnan(Lmax))~=OCSHOOTCONT.TIMEMESH.num % if arcid cannot be determined for every time point return empty
            ctrl_loc=[];
            return
        end

        % determine control values depending on the arcid
        ctrl_loc=OCGRADSOL.optimalcontrol(t,depvar,par,arcid);
    case 'posterior'
        %
end

function adm=isadmissible(t,depvar,ctrl_locnew,par,arcid)
global OCSHOOTCONT OCGRADSOL

adm=OCGRADSOL.admissible(t,depvar,ctrl_locnew,par,arcid);

