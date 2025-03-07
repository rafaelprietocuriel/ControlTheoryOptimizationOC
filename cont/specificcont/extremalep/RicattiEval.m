% Ricatti evaluation file
% =======================

function result=RicattiEval(xhat,Y,modelpar,arc,unstable_flag,Y)

global OCMATAE
J=OCMATAE.niteration;
A = hetjac(xhat,p,J);
result = [];

if unstable_flag
  if  ~isempty(Y)
    Q0U = OCMATAE.Q0;
    % Riccati blocks from unstable eigenspace
    [U11, U12, UE21, U22]=RicattiCoeff(Q0U,A,OCMATAE.nu);
    tmp = (U22*Y - Y*U11 + UE21 - Y*U12*Y)';
    for i=1:OCMATAE.ns
        result(end+1:end+OCMATAE.nu,1) = tmp(:,i);
    end
  end
else
  if  ~isempty(Y)
    Q1S = OCMATAE.Q1;
    % Riccati blocks from stable eigenspace
    [S11, S12, SE21, S22] = Het_RicattiCoeff(Q1S,A,OCMATAE.ns);
    tmp = (S22*Y - Y*S11 + SE21 - Y*S12*Y)';
    for i=1:OCMATAE.nu
        result(end+1:end+OCMATAE.ns,1) = tmp(:,i);
    end
  end
end
