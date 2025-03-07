function disp(ocgLim)

for ii=1:multiplicity(ocgLim)
    if multiplicity(ocgLim)>1
        disp(['Limitset ' num2str(ii)])
    end
    equilibflag=isequilibrium(ocgLim(ii));
    eigval=eig(ocgLim(ii));
    if isequal(get(0,'FormatSpacing'),'compact')
        disp([inputname(1) ' =']);
        if ~isempty(ocgLim(ii))
            disp('ocmatclass: ocgradlimitset')
            disp(['modelname: ' modelname(ocgLim(ii))])
            if equilibflag
                disp('Equilibrium:')
            else
                disp('Periodic Solution:')
            end
            disp(state(ocgLim(ii)));
            disp(costate(ocgLim(ii)));
        else
            disp('empty ocgradlimitset');
        end
    else
        disp(' ')
        disp([inputname(1) ' =']);
        disp(' ')
        if ~isempty(ocgLim(ii))
            disp('ocmatclass: ocgradlimitset')
            disp(' ');
            disp(['modelname: ' modelname(ocgLim(ii))])
            if equilibflag
                disp(' ');
                disp('Equilibrium:')
            else
                disp(' ');
                disp('Period:')
                disp(period(ocgLim(ii)))
                disp(' ');
                disp('Periodic Solution:')
            end
            disp(' ');
            disp('State:')
            disp(state(ocgLim(ii)));
            disp(' ');
            disp('Costate:')
            disp(costate(ocgLim(ii)));
            disp(' ');
            disp('Eigenvalues:')
            if isempty(eigval)
                disp('  []')
            else
                disp(eigval);
            end
            disp(' ');
        else
            disp('empty ocgradlimitset');
            disp(' ');
        end
    end
end