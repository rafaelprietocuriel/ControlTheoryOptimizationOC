function [arcpart4yidx,arcpart4fzidx,arcpart4fFintidx,arcpart4zFintidx,arcpart4tFintidx]=calcresarcpartidx_dae(Nm1)
% 
% for the computation the solution components are reordered such that the
% first order components are at the beginnig and thereafter the zero order
% components. The sol.y components are in the initial ordering.

global OCMATCONT
persistent N arcpart4ycoefficientidx arcpart4firstorderzcoefficientidx  ...
    arcpart4firstorderFinterioridx arcpart4zeroorderFinterioridx arcpart4transitionFidx %arcpart4zeroorderidx

if isempty(N) || all(N~=Nm1)
    N=Nm1;
    if length(N)==1
        arcpart4ycoefficientidx=[1;N];
        arcpart4firstorderzcoefficientidx=[1;OCMATCONT.CollocationNumber*N];
        %arcpart4zeroorderidx=[1;OCMATCONT.zeroordercomponentnumber*N];
        arcpart4firstorderFinterioridx=[1;OCMATCONT.firstordercomponentnumber*OCMATCONT.CollocationNumber*N];
        arcpart4zeroorderFinterioridx=[1;OCMATCONT.zeroordercomponentnumber*OCMATCONT.CollocationNumber*N];
        arcpart4transitionFidx=[1;OCMATCONT.firstordercomponentnumber*(N-1)];

    else
        arcpart4ycoefficientidx=cumsum(N);
        arcpart4ycoefficientidx=[1 arcpart4ycoefficientidx(1:OCMATCONT.arcnumber-1)+1;arcpart4ycoefficientidx];

        arcpart4firstorderzcoefficientidx=cumsum(OCMATCONT.CollocationNumber*N);
        arcpart4firstorderzcoefficientidx=[1 arcpart4firstorderzcoefficientidx(1:OCMATCONT.arcnumber-1)+1;arcpart4firstorderzcoefficientidx];

        %arcpart4zeroorderidx=cumsum(OCMATCONT.zeroordercomponentnumber*N);
        %arcpart4zeroorderidx=[1 arcpart4zeroorderidx(1:OCMATCONT.arcnumber-1)+1;arcpart4zeroorderidx];

        arcpart4firstorderFinterioridx=OCMATCONT.firstordercomponentnumber*OCMATCONT.CollocationNumber*cumsum(OCMATCONT.meshNumber-1);
        arcpart4firstorderFinterioridx=[1 arcpart4firstorderFinterioridx(1:end-1)+1;arcpart4firstorderFinterioridx];
        
        arcpart4zeroorderFinterioridx=OCMATCONT.zeroordercomponentnumber*OCMATCONT.CollocationNumber*cumsum(OCMATCONT.meshNumber-1);
        arcpart4zeroorderFinterioridx=[1 arcpart4zeroorderFinterioridx(1:OCMATCONT.arcnumber-1)+1;arcpart4zeroorderFinterioridx];

        arcpart4transitionFidx=OCMATCONT.firstordercomponentnumber*cumsum(N-1);
        arcpart4transitionFidx=[1 arcpart4transitionFidx(1:OCMATCONT.arcnumber-1)+1;arcpart4transitionFidx];
    end
end
arcpart4yidx=arcpart4ycoefficientidx;
arcpart4fzidx=arcpart4firstorderzcoefficientidx;
%arcpart4zidx=arcpart4zeroorderidx;
arcpart4fFintidx=arcpart4firstorderFinterioridx;
arcpart4zFintidx=arcpart4zeroorderFinterioridx;
arcpart4tFintidx=arcpart4transitionFidx;