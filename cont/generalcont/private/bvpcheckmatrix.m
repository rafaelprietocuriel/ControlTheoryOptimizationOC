function isSingular = bvpcheckmatrix(A,L,U,P,Q,R)
%BVPCHECKMATRIX  Helper function to check A (=R*P'*L*U*Q') for 'near singularity'
%
%   See also BVP4C, BVP5C.

%   Copyright 2007 The MathWorks, Inc.
%   $Revision: 1.1.6.2 $  $Date: 2007/06/14 05:09:00 $

% Adjust the warnings.
global OCBVP
[lastmsg,lastid] = lastwarn('');
warning(OCBVP.warnoff);

Ainv_norm = normest1(@condaux);

warning(OCBVP.warnstat);

[msg,msgid] = lastwarn;
if isempty(msg)
    lastwarn(lastmsg,lastid);
elseif any(strcmp(msgid,OCBVP.warnoffId))
    isSingular = true;
    return;
end

%isSingular = (Ainv_norm*norm(A,inf)*eps > 1e-5);
isSingular = (Ainv_norm*norm(A,inf)*eps > 1e-2);

% --------------------------------------------------------
% Nested function, with access to L,U,P,Q,R
% --------------------------------------------------------

    function f = condaux(flag,X)
        %CONDAUX  Auxiliary function for estimation of condition.
        switch flag
            case 'dim'
                f = max(size(L));
            case 'real'
                f = 1;
            case 'notransp'
                f = Q * (U \ (L \ (P * (R \ X))));
            case 'transp'
                f = R' \ (P' * (L' \ (U' \ (Q' * X))));
        end
    end  % condaux

% --------------------------------------------------------

end  % bvpsingularmatrix
