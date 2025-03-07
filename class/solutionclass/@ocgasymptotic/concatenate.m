function ocgAsym=concatenate(varargin)

if nargin==2
    ocTrj1=ocgtrajectory(varargin{1});
    ocTrj2=ocgtrajectory(varargin{2});
    ocTrj=concatenate(ocTrj1,ocTrj2);
    ocgAsym=ocgasymptotic(ocTrj,limitset(varargin{2}));
else
    ocgAsym=concatenate(varargin{1},concatenate(varargin{2:end}));
end