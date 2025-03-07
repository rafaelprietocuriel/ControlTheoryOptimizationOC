function initocmatcontdomainvariable(solver,domaindata,discretizationdata)
global OCMATCONT

for ii=1:numel(domaindata)
    OCMATCONT.DOMAINDDATA(ii).numode=domaindata(ii).odedim;
    OCMATCONT.DOMAINDDATA(ii).numae=domaindata(ii).aedim;
    OCMATCONT.DOMAINDDATA(ii).numeq=numel(domaindata(ii).daeorder);%number of equations
    OCMATCONT.DOMAINDDATA(ii).eqcoord=1:OCMATCONT.DOMAINDDATA(ii).numeq;
    OCMATCONT.DOMAINDDATA(ii).odecoord=1:domaindata(ii).odedim;
    OCMATCONT.DOMAINDDATA(ii).aecoord=domaindata(ii).odedim+(1:domaindata(ii).aedim);
end

switch solver
    case 'sbvpoc'
        OCMATCONT.solver='sbvpoc';
        for ii=1:numel(domaindata)
            [OCMATCONT.DOMAINDDATA(ii).nodes,OCMATCONT.DOMAINDDATA(ii).psi, ...
                OCMATCONT.DOMAINDDATA(ii).psival, ...
                OCMATCONT.DOMAINDDATA(ii).weights, ...
                OCMATCONT.DOMAINDDATA(ii).psi0, ...
                OCMATCONT.DOMAINDDATA(ii).psival0, ...
                OCMATCONT.DOMAINDDATA(ii).psivalequi, ...
                OCMATCONT.DOMAINDDATA(ii).psivalequi0]=discretizationbasisdata(max([domaindata(ii).daeorder]),discretizationdata(ii).collocationmethod,discretizationdata(ii).numcollocationpoints);
            %OCMATCONT.DOMAINDDATA(ii).psival2=permute(OCMATCONT.DOMAINDDATA(ii).psival,[3 1 2]);
            OCMATCONT.DOMAINDDATA(ii).method=discretizationdata(ii).collocationmethod;
            OCMATCONT.DOMAINDDATA(ii).order=domaindata(ii).daeorder;
            OCMATCONT.DOMAINDDATA(ii).maxorder=max(domaindata(ii).daeorder);
            OCMATCONT.DOMAINDDATA(ii).minorder=min(domaindata(ii).daeorder);
            OCMATCONT.DOMAINDDATA(ii).numode=domaindata(ii).odedim;
            OCMATCONT.DOMAINDDATA(ii).numae=domaindata(ii).aedim;
            OCMATCONT.DOMAINDDATA(ii).numeq=numel(domaindata(ii).daeorder);%number of equations
            OCMATCONT.DOMAINDDATA(ii).eqcoord=1:OCMATCONT.DOMAINDDATA(ii).numeq;
            OCMATCONT.DOMAINDDATA(ii).odecoord=1:domaindata(ii).odedim;
            % it is assumed that algebraic equations appear at the end of
            % the DAEs
            OCMATCONT.DOMAINDDATA(ii).aecoord=domaindata(ii).odedim+(1:domaindata(ii).aedim);
            OCMATCONT.DOMAINDDATA(ii).numcols=discretizationdata(ii).numcollocationpoints;
            OCMATCONT.DOMAINDDATA(ii).numcolsp1=OCMATCONT.DOMAINDDATA(ii).numcols+1;
            OCMATCONT.DOMAINDDATA(ii).numcolscoord=1:OCMATCONT.DOMAINDDATA(ii).numcols;
            OCMATCONT.DOMAINDDATA(ii).numeqnumcols=OCMATCONT.DOMAINDDATA(ii).numeq*OCMATCONT.DOMAINDDATA(ii).numcols;
            OCMATCONT.DOMAINDDATA(ii).numeqnumcolscoord=1:OCMATCONT.DOMAINDDATA(ii).numeqnumcols;
        end
        %OCMATCONT.factorial=[1;1;2;6;24;120;720;5040;40320;362880;3628800;39916800;479001600;6227020800];
        warnstat(1) = warning('query','MATLAB:singularMatrix');
        warnstat(2) = warning('query','MATLAB:nearlySingularMatrix');
        warnoff = warnstat;
        warnoff(1).state = 'off';
        warnoff(2).state = 'off';
        OCMATCONT.WARNING.warnstat=warnstat;
        OCMATCONT.WARNING.warnoff=warnoff;
    case {'bvp4c','bvp6c'}
        OCMATCONT.solver=solver;
        
        % bvpsolver specific variables
        warnstat(1) = warning('query','MATLAB:singularMatrix');
        warnstat(2) = warning('query','MATLAB:nearlySingularMatrix');
        warnoff = warnstat;
        warnoff(1).state = 'off';
        warnoff(2).state = 'off';
        OCMATCONT.BVPSOLVERVARIABLE.(solver).warnstat=warnstat;
        OCMATCONT.BVPSOLVERVARIABLE.(solver).warnoff=warnoff;
        % for bvp4c and bvp6c the number of ODEs must not change for
        % different arcs
        OCMATCONT.BVPSOLVERVARIABLE.(solver).numode=OCMATCONT.DOMAINDDATA(1).odedim;
        OCMATCONT.BVPSOLVERVARIABLE.(solver).numae=OCMATCONT.DOMAINDDATA(1).aedim;
        OCMATCONT.BVPSOLVERVARIABLE.(solver).numeq=OCMATCONT.DOMAINDDATA(1).numeq;
        OCMATCONT.BVPSOLVERVARIABLE.(solver).eqcoord=OCMATCONT.DOMAINDDATA(ii).eqcoord;
        OCMATCONT.BVPSOLVERVARIABLE.(solver).odecoord=OCMATCONT.DOMAINDDATA(ii).odecoord;
    case 'bvp5c'
        
        OCMATCONT.solver=solver;
        for ii=1:numel(domaindata)
            OCMATCONT.DOMAINDDATA(ii).numode=domaindata(ii).odedim;
            OCMATCONT.DOMAINDDATA(ii).numae=domaindata(ii).aedim;
            OCMATCONT.DOMAINDDATA(ii).numeq=numel(domaindata(ii).daeorder);%number of equations
            OCMATCONT.DOMAINDDATA(ii).eqcoord=1:OCMATCONT.DOMAINDDATA(ii).numeq;
            OCMATCONT.DOMAINDDATA(ii).odecoord=1:domaindata(ii).odedim;
        end
        
        % bvpsolver specific variables
        warnstat(1) = warning('query','MATLAB:singularMatrix');
        warnstat(2) = warning('query','MATLAB:nearlySingularMatrix');
        warnoff = warnstat;
        warnoff(1).state = 'off';
        warnoff(2).state = 'off';
        OCMATCONT.BVPSOLVERVARIABLE.(solver).warnstat=warnstat;
        OCMATCONT.BVPSOLVERVARIABLE.(solver).warnoff=warnoff;
        % for bvp4c and bvp6c the number of ODEs must not change for
        % different arcs
        OCMATCONT.BVPSOLVERVARIABLE.(solver).numode=OCMATCONT.DOMAINDDATA(1).odedim;
        OCMATCONT.BVPSOLVERVARIABLE.(solver).numae=OCMATCONT.DOMAINDDATA(1).aedim;
        OCMATCONT.BVPSOLVERVARIABLE.(solver).numeq=OCMATCONT.DOMAINDDATA(1).numeq;
        OCMATCONT.BVPSOLVERVARIABLE.(solver).eqcoord=OCMATCONT.DOMAINDDATA(ii).eqcoord;
        OCMATCONT.BVPSOLVERVARIABLE.(solver).odecoord=OCMATCONT.DOMAINDDATA(ii).odecoord;
end
