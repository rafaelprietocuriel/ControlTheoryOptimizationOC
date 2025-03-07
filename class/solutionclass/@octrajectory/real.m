function realocTrj=real(ocTrj)
%
% REAL returns the real valued counterpart.
%
% RDYNPRIM = REAL(DYNPRIM) returns a dynprimitive object RDYNPRIM, with only 
% the real parts, as well for the dynVar field as the linearization field,
% of the dynprimitive DYNPRIM.

realocTrj=ocTrj;
realocTrj.y=real(realocTrj.y);
realocTrj.arcinterval=real(realocTrj.arcinterval);
realocTrj.modelparameter=real(realocTrj.modelparameter);
realocTrj.timehorizon=real(realocTrj.timehorizon);
realocTrj.linearization=real(realocTrj.linearization);
fn=fieldnames(realocTrj.solverinfo);
for ii=1:length(fn)
    if isnumeric(realocTrj.solverinfo.(fn{ii}))
        realocTrj.solverinfo.(fn{ii})=real(realocTrj.solverinfo.(fn{ii}));
    end
end