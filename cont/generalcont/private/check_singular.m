function isSingular = check_singular(A,L,U,P)
%CHECK_SINGULAR  Check A (L*U*P) for 'singularity'; mute certain warnings.

global OCBVP OCMATCONT
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

isSingular = (Ainv_norm*norm(A,inf)*eps > OCMATCONT.OPTIONS.singularthreshold);

%--------------------------------------------------------------------------

    function f = condaux(flag,X)
        %CONDAUX  Auxiliary function for estimation of condition.

        switch flag
            case 'dim'
                f = max(size(L));
            case 'real'
                f = 1;
            case 'notransp'
                f = U \ (L \ (P * X));
            case 'transp'
                f = P * (L' \ (U' \ X));
        end
    end
end

%--------------------------------------------------------------------------



