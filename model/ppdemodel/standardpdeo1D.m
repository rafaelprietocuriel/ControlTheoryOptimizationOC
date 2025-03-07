classdef standardpdeo1D < pde 
% standardpdeo1D:  standard 1D OOPDE object for ocmat (classdef) 
methods(Access = public)
  function o=standardpdeo1D(varargin) % constructor 
      o.grid=grid1D; 
      
      if nargin==2
          lb=varargin{1};
          la=-lb;
          h=varargin{2};
      elseif nargin==3
          la=varargin{1};
          lb=varargin{2};
          h=varargin{3};
      end
      
      o.grid.interval([la,lb],h); 
      % test functions are given by the Lagrange functions
      o.fem=lagrange11D;
  end
end  
methods(Access = protected) % only here since required by pde-class
    function r=df(~,~,~)
        r=0; 
    end
end
end
