%-----------------------------------------------------------------
function [dFdy,nfevals,nfcalls]=numjaccsd(fun,Fargs,nF,options)
% NUMJACCSD    Complex Step Jacobian
% based on 'jacobiancsd' by Yi Cao at Cranfield University, 02/01/2008 and
% uses the structure of odenumjac.
%
% J = NUMJACCSD(F,FARGS,NF,OPTIONS) returns the numerical (NF x N) Jacobian
% matrix of a NF-vector function, F(FARGS{:}) at the reference point, 
% Y=FARGS{OPTIONS.DIFFVAR} (N-vector). 
%
% The structure OPTIONS must have the following fields: DIFFVAR, VECTVARS.
% The field OPTIONS.DIFFVAR is the index of the  differentiation variable,
% Y = FARGS{DIFFVAR}. For a function F(t,x), set DIFFVAR to 1 to compute
% DF/Dt, or to 2 to compute DF/Dx. Set OPTIONS.VECTVAR to the indices of
% vectorized arguments: VECTVAR = [2] indicates that  F(t,[x1 y2 ...])
% returns [F(t,x1) F(t,x2) ...], while VECTVAR = [1,2] indicates that F([t1
% t2 ...],[x1 x2 ...]) returns [F(t1,x1) F(t2,x2) ...].
%   
% [DFDY,NFEVALS,NFCALLS] = ODENUMJAC(...) returns the number of values
% F(FARGS{:}) computed while forming dFdY (NFEVALS) and the number of calls
% to the function F (NFCALLS). If F is not vectorized, NFCALLS equals
% NFEVALS.  
% Dieter Grass

% Options
diffvar = options.diffvar; 
vectvar = options.vectvars;

% The differentiation variable.
y  = Fargs{diffvar};
classY = class(y);

ny = length(y);

dFdy=zeros(nF,ny);                   % allocate memory for the Jacobian matrix
del = (y + ny*eps(classY)*i) - y;
h=imag(del);
ydel = y(:,ones(1,ny)) + diag(del);
if isempty(vectvar)
    for ii=1:ny                      % loop for each independent variable
        dFdy(:,ii)=imag(fun(Fargs{1:diffvar-1},ydel(:,ii),Fargs{diffvar+1:end}))/h(ii);     % complex step differentiation
    end
    nfcalls = ny;                       % stats
else
    Fargs_expanded = Fargs;
    Fargs_expanded{diffvar} = ydel;
    vectvar = setdiff(vectvar,diffvar);
    for ii=1:length(vectvar)
      Fargs_expanded{vectvar(ii)} = repmat(Fargs{vectvar(ii)},1,ny);
    end
    dFdy=imag(fun(Fargs_expanded{:}))./h(:,ones(1,ny)); 
    nfcalls = 1;                       % stats
end
nfevals = ny;                         % stats (at least one per loop)
