function ocComp=occomposite(varargin)
%
% OCCOMPOSITE occomposite constructor
%
% OCCOMPOSITE
%

%

switch nargin
    case 0
        ocComp.path{1}=[];
        ocComp.order=[];
        ocComp.userinfo=struct;
        ocComp.label='';
        ocComp=class(ocComp,'occomposite');

    otherwise
        ocComp=[];
        counter=0;
        for ii=1:nargin
            if  ~isempty(varargin{ii})&& (isocasymptotic(varargin{ii}) || isoctrajectory(varargin{ii}) || ismmultipath(varargin{ii}))
                counter=counter+1;
                ocComp.path{counter}=varargin{ii};
            elseif iscell(varargin{1})
                ocComp=occomposite(varargin{1}{:});
                return
            end
        end

        if isempty(ocComp)
            ocComp=occomposite;
            return
        end
        ocComp.order=counter;
        ocComp.userinfo=struct;
        ocComp.label='';
        ocComp=class(ocComp,'occomposite');
end
