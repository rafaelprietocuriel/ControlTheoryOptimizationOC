function [varargout] = gradmap_old(u,v,w,FGRAD,options,varargin)
%FMINGRADMAPOPT Finds a minimum of a function of several variables constrained on a simple set.
%   FMINCON solves problems of the form:
%       min F(X)  subject to:  LB <= X <= UB  (Box constraints)           
%        X                     |X| <= r       (Euclidean Ball constraints)
%                                
%   by Gradient Mapping Optimal Method for convex function on simple sets
%   this is the constant step scheme

% check input parameters
if nargin < 5
  options = [];
  if nargin < 4
      error('gradmap requires at least four input arguments');
end,end

% read options
Proj = 0;% Set control bounds to NaN if they are not cpecified
if ~isempty(u) 
  if ~isfield(options,'LowerBoundU'),options.LowerBoundU = NaN;
  elseif isempty(options.LowerBoundU),error('LowerBoundU should be nonemty array, scalar or NaN')
  else Proj = 1; 
  end
  if ~isfield(options,'UpperBoundU'),options.UpperBoundU = NaN;
  elseif isempty(options.UpperBoundU),error('UpperBoundU should be nonemty array, scalar or NaN')
  else Proj = 1; 
  end
end
if ~isempty(v)
  if ~isfield(options,'LowerBoundV'),options.LowerBoundV = NaN;
  elseif isempty(options.LowerBoundV),error('LowerBoundV should be nonemty array, scalar or NaN')
  else Proj = 1; 
  end
  if ~isfield(options,'UpperBoundV'),options.UpperBoundV = NaN;
  elseif isempty(options.UpperBoundV),error('UpperBoundV should be nonemty array, scalar or NaN')
  else Proj = 1; 
  end
end
if ~isempty(w) 
  if ~isfield(options,'LowerBoundW'),options.LowerBoundW = NaN;
  elseif isempty(options.LowerBoundW),error('LowerBoundW should be nonemty array, scalar or NaN')
  else Proj = 1; 
  end
  if ~isfield(options,'UpperBoundW'),options.UpperBoundW = NaN;
  elseif isempty(options.UpperBoundW),error('UpperBoundW should be nonemty array, scalar or NaN')
  else Proj = 1; 
  end
end
if isfield(options,'Projection')
  if Proj == 1, warning('gradmap:ControlBounds','Control Bounds might be inefective since Projection is specified!')
  end
  if isa(options.Projection, 'function_handle'), Proj = 4;
  elseif ischar(options.Projection), Proj = 3;
  else error('Projection: ''PositiveOrthant'',''UnitBox'', or function_handle')
  end
elseif Proj == 1, options.Projection = @ValueBoundBox; %@UserDefinedBox;
end

DisplayGrad = 0;
DisplayVal = 0;
if isfield(options,'Display')
    if ischar(options.Display) 
        if strcmp(options.Display,'iter') 
            DisplayGrad = 1;
            DisplayVal = 1;
        elseif strcmp(options.Display,'iterVal') 
            DisplayVal = 1;
        elseif strcmp(options.Display,'iterGrad') 
            DisplayGrad = 1;
        end 
        % print caption
        if DisplayGrad == 1 || DisplayVal == 1
          fprintf('Iter ');
          if DisplayVal == 1, fprintf('\t ObjVal '); end
          if DisplayGrad == 1
            if Proj % if there is a projection we display projected gradient
              if ~isempty(u),fprintf('\t Max(Proj(dU)) ');end
              if ~isempty(v),fprintf('\t Max(Proj(dV)) ');end
              if ~isempty(w),fprintf('\t Max(Proj(dW)) ');end
            else
              if ~isempty(u),fprintf('\t Max(dU) ');end
              if ~isempty(v),fprintf('\t Max(dV) ');end
              if ~isempty(w),fprintf('\t Max(dW) ');end
            end
          end
          fprintf('\n');
        end
    end
end

LineSrch = 0;
if isfield(options,'LineSrch')
    if ischar(options.LineSrch) 
        if strcmp(options.LineSrch,'GradChord'), LineSrch = 1;    
        elseif strcmp(options.LineSrch,'ValCubic') 
            LineSrch = -1;
            gmopt = struct('CalcVal',1, 'CalcGrad',0, 'DisplayVal',0);
        elseif strcmp(options.LineSrch,'NONE'), LineSrch = 0;
        else
          error('There is no such LineSrch option! Only GradChord, ValCubic or NONE');
        end 
    end
end

MaxIter = 100;% Maximal iterations.
if isfield(options,'MaxIter')
    if isreal(options.MaxIter), MaxIter = options.MaxIter; 
    end
end

L = 100;% Lipschitz constatnt. If L < 0 then it is a maximization
if isfield(options,'LipschitzConst')
    if isreal(options.LipschitzConst), L = options.LipschitzConst; 
    end
end

% Minimize or Maximize. if it is not set the sign of L determines that
% If L < 0 then it is a maximization
if isfield(options,'Optimization')
    if istrcmp(options.Optimization,'max'), L = -abs(L); 
    elseif istrcmp(options.Optimization,'min'), L = abs(L); 
    end
end

ConditionNumber = 1; % Condition number is the ratio of maximal second derivative to minimal second derivative
if isfield(options,'ConditionNumber')
    if isreal(options.ConditionNumber), ConditionNumber = options.ConditionNumber; 
    end
end
if ConditionNumber < 1, error('Condition Number of the function cannot be less than 1'); end

b = sqrt(ConditionNumber); % parameter for constant step scheeme
b = (b-1)/(b+1);

% n = 1; % Number of arguments
% if isfield(options,'ArgNumber')
%     if isreal(options.ArgNumber), n = options.ArgNumber; 
% end,end
% if n < 1, error('Number of arguments the function must be grater than 0'); end

% make an initial projection
if Proj
[u,v,w] = feval(options.Projection,u,v,w,options,varargin{:});
end

% auxiliary vector Y which can be out of the set
yu = u;
yv = v;
yw = w;
% auxiliary vector which keeps the previous position
u0 = u;
v0 = v;
w0 = w;
% auxiliary vector which keeps the gradient
du = u;
dv = v;
dw = w;
%
impr = 1;
Iter = 0;
Step = 1/L;
while Iter < MaxIter %% && abs(Step) >= 1/abs(L)/4
    Iter = Iter + 1;
    %
    if LineSrch < 0 || DisplayVal
      [yu,yv,yw,y0,du,dv,dw] = feval(FGRAD,yu,yv,yw,...
         struct('CalcVal',1, 'CalcGrad',1, 'DisplayVal',0),varargin{:});
    else
      [yu,yv,yw,y0,du,dv,dw] = feval(FGRAD,yu,yv,yw,[],varargin{:}); % y0 = NaN
    end
    if DisplayGrad == 1 || DisplayVal == 1
      fprintf('%g ',Iter);
      if DisplayVal == 1, fprintf('\t %16.8f',y0); end
      if DisplayGrad == 1
        if Proj % if there is a projection we display projected gradient
          [u,v,w] = feval(options.Projection,yu-Step*du,yv-Step*dv,yw-Step*dw,options,varargin{:});
          if ~isempty(du),fprintf('\t %16.8f',max(abs(du(:).*(u(:)~=yu(:)))));end
          if ~isempty(dv),fprintf('\t %16.8f',max(abs(dv(:).*(v(:)~=yv(:)))));end
          if ~isempty(dw),fprintf('\t %16.8f',max(abs(dw(:).*(w(:)~=yw(:)))));end
        else
          if ~isempty(du),fprintf('\t %16.8f',max(abs(du(:))));end
          if ~isempty(dv),fprintf('\t %16.8f',max(abs(dv(:))));end
          if ~isempty(dw),fprintf('\t %16.8f',max(abs(dw(:))));end
        end
      end
      fprintf('\n');
    end
%     if Iter >= 536
%         fprintf('\t %16.8f',y0);
%     end
    % get new variables
    u = yu - Step * du;
    v = yv - Step * dv;
    w = yw - Step * dw;
    % adjust the step using gradients or value function
    if LineSrch > 0 % adjust the step using gradients by chord method
      if Proj, [u,v,w] = feval(options.Projection,u,v,w,options,varargin{:});end
      % calculate the squared norm of the gradient
      y1 = [du(:);dv(:);dw(:)]' * [du(:);dv(:);dw(:)];
      y2 = 1;
      s1 = 0;
      s2 = Step;
      s = Step;
      impr = 0;
      while Iter < MaxIter % && abs(Step) >= 1/abs(L)/4 && (Prod2 < 0 || Norm1 < Prod2)
        Iter = Iter + 1;
        % we write gradient into the variables u,v,w to economize memory
        [u,v,w,y0,du2,dv2,dw2] = feval(FGRAD,u,v,w,[],varargin{:});
        % calculate scalar product of two gradients
        y = [du(:);dv(:);dw(:)]' * [du2(:);dv2(:);dw2(:)];
        if abs(y) < 1/abs(L)
            Step = s;
            impr = 1;
            break              % if y is of different sign with y1 or y2
        elseif y < 0 && (y1 > 0 || y2 > 0) || y > 0 && (y1 < 0 || y2 < 0)
            if y < 0 && y1 > 0, y2 = y; s2 = s; 
            elseif y < 0 && y2 > 0, y1 = y; s1 = s; 
            elseif y > 0 && y1 < 0, y2 = y; s2 = s; 
            elseif y > 0 && y2 < 0, y1 = y; s1 = s; 
            end % chord aproximation
            s = s1 - y1 * (s2 - s1) / (y2 - y1); % chord aproximation     
        elseif y1 <= y % if not convexity (concavity if L<0)
            Step = Step * 1.5; 
            impr = 1; break
        else impr = 1; break
        end        
        u = yu - s * du;
        v = yv - s * dv;
        w = yw - s * dw;
        if Proj, [u,v,w] = feval(options.Projection,u,v,w,options,varargin{:});end
      end
      if impr ~= 1
        u = u0;
        v = v0;
        w = w0;
        if DisplayVal == 1, fprintf('%g \t No improvement \n',Iter); end
        break;
      end
      %if Proj, [u,v,w] = feval(options.Projection,u,v,w,options,varargin{:});end
    elseif LineSrch < 0 % adjust the step using value function
      impr = 0;
      while Iter < MaxIter && impr == 0
        Iter = Iter + 1;
        x1 = Step;
        if Proj, [u,v,w] = feval(options.Projection,u,v,w,options,varargin{:});end
        [u,v,w,y1] = feval(FGRAD,u,v,w,gmopt,varargin{:});
        if ((y1 < y0) && (L > 0)) || ((y1 > y0) && (L < 0))
          OFVal = y1;
          impr = 1;
          Step = 1.5 * Step;
          x2 = Step;
          Ubest = u;
          Vbest = v;
          Wbest = w;
          u = yu - Step * du;
          v = yv - Step * dv;
          w = yw - Step * dw;
          if Proj, [u,v,w] = feval(options.Projection,u,v,w,options,varargin{:});end
          [u,v,w,y2] = feval(FGRAD,u,v,w,gmopt,varargin{:});
          Iter = Iter + 1;
          if ((y2 < y1) && (L > 0)) || ((y2 > y1) && (L < 0))
            if (x1*y2 + x2*y0 - (x1+x2)*y1 >= 0)%does not depend on sign(L)
              OFVal = y2;
            else
              Ubest = u;
              Vbest = v;
              Wbest = w;
              s = ((x2^2 - x1^2)*y0 + y2*x1^2 - y1*x2^2)/(2*(y2*x1 - y1*x2 + y0*(x2-x1)));  %[x1,x2,s]
              u = yu - s * du;
              v = yv - s * dv;
              w = yw - s * dw;
              if Proj, [u,v,w] = feval(options.Projection,u,v,w,options,varargin{:});end
              [u,v,w,y3] = feval(FGRAD,u,v,w,gmopt,varargin{:});
              Iter = Iter + 1;
              if ((y3 < y2) && (L > 0)) || ((y3 > y2) && (L < 0)) && ~isinf(s) %[1,y0,y1,y2,OFValnew]    
                Step = s;
                OFVal = y3;
              else
                u = Ubest;
                v = Vbest;
                w = Wbest;
                OFVal = y2;
              end
            end   
          else
            if (y2*x1 - y1*x2 + y0*(x2-x1) <= 0)%does not depend on sign(L)
              u = Ubest;
              v = Vbest;
              w = Wbest;
              OFVal = y1;
            else
              s = ((x2^2 - x1^2)*y0 + y2*x1^2 - y1*x2^2)/(2*(y2*x1 - y1*x2 + y0*(x2-x1)));  %[x1,x2,s] 
              u = yu - s * du;
              v = yv - s * dv;
              w = yw - s * dw;
              if Proj, [u,v,w] = feval(options.Projection,u,v,w,options,varargin{:});end
              [u,v,w,y3] = feval(FGRAD,u,v,w,gmopt,varargin{:});
              Iter = Iter + 1;
              if ((y3 < y1) && (L > 0)) || ((y3 > y1) && (L < 0)) % !!!! %[2,y0,y1,y2,OFValnew]
                OFVal = y3;
                Step = s;
              else
                u = Ubest;
                v = Vbest;
                w = Wbest;
                OFVal = y1;
                Step = x1;
              end
            end
          end   
        else   %(y1 < y2)
          Step = 0.5 * Step;
          x2 = Step;
          u = yu - Step * du;
          v = yv - Step * dv;
          w = yw - Step * dw;
          if Proj, [u,v,w] = feval(options.Projection,u,v,w,options,varargin{:});end
          [u,v,w,y2] = feval(FGRAD,u,v,w,gmopt,varargin{:});
          Iter = Iter + 1;
          if ((y2 < y0) && (L > 0)) || ((y2 > y0) && (L < 0))
            impr = 1;
            if (y2*x1 - y1*x2 + y0*(x2-x1) >= 0)%does not depend on sign(L)
              OFVal = y2;
            else
              Ubest = u;
              Vbest = v;
              Wbest = w;
              s = ((x2^2 - x1^2)*y0 + y2*x1^2 - y1*x2^2)/(2*(y2*x1 - y1*x2 + y0*(x2-x1)));   %[x1,x2,s]
              u = yu - s * du;
              v = yv - s * dv;
              w = yw - s * dw;
              if Proj, [u,v,w] = feval(options.Projection,u,v,w,options,varargin{:});end
              [u,v,w,y3] = feval(FGRAD,u,v,w,gmopt,varargin{:});
              Iter = Iter + 1;
              if ((y3 < y2) && (L > 0)) || ((y3 > y2) && (L < 0))  %  [3,y0,y1,y2,OFValnew]
                OFVal = y3;
                Step = s;
              else
                u = Ubest;
                v = Vbest;
                w = Wbest;
                OFVal = y2;
              end
            end
          elseif ((y2 < 0.5 * (y0 + y1)) && (L > 0)) || ((y2 > 0.5 * (y0 + y1)) && (L < 0)) 
            s = ((x2^2 - x1^2)*y0 + y2*x1^2 - y1*x2^2)/(2*(y2*x1 - y1*x2 + y0*(x2-x1)));  %[x1,x2,s] check y2*x1 - y1*x2 + y0*(x2-x1) >=0 
            if ((s > 0) && (L > 0) || (s < 0) && (L < 0)) && ~isinf(s) %!!! % check for s equal inf if y2*x1 - y1*x2 + y0*(x2-x1) == 0
              u = yu - s * du;
              v = yv - s * dv;
              w = yw - s * dw;
              if Proj, [u,v,w] = feval(options.Projection,u,v,w,options,varargin{:});end
              [u,v,w,y3] = feval(FGRAD,u,v,w,gmopt,varargin{:});
              Iter = Iter + 1;
              if ((y3 <= y0) && (L > 0)) || ((y3 >= y0) && (L < 0))   % [4,y0,y1,y2,OFValnew]
                impr = 1;
                OFVal = y3;
                Step = s;
              end
            end
          end       
        end
        if DisplayVal == 1
          if impr == 1 
           %fprintf('%g Success: ObjVal =  %16.8f \n', [Iter,  OFVal] );
          else
           if Iter < MaxIter, fprintf('%g \t No improvement \n',Iter)
           end
          end
        end
        if impr ~= 1
          Step = 0.5 * Step;
          if abs(Step) < 1/abs(L)/4, break; end
        end
      end
    elseif Proj, [u,v,w] = feval(options.Projection,u,v,w,options,varargin{:});
    end
    if impr ~= 1
        u = u0;
        v = v0;
        w = w0;
        break;
    end
    % calculate new Y
    yu = u + b * (u - u0);
    yv = v + b * (v - v0);
    yw = w + b * (w - w0);
    % for convex (concave) function this projection is unnesessery
    % if Proj, [yu,yv,yw] = feval(options.Projection,yu,yv,yw,options,varargin{:});end
    % remember the previous position
    u0 = u;
    v0 = v;
    w0 = w;
end
[varargout{1:nargout(func2str(FGRAD))}] = feval(FGRAD,u,v,w,...
         struct('CalcVal',1, 'CalcGrad',1, 'DisplayVal',0),varargin{:});
end
