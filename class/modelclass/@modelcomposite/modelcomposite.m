function mmObj=modelcomposite(varargin)
%
% MODELCOMPOSITE initializes a class combining multiple models. This object 
% consists of the two fields MODEL and RESULT. 

switch nargin
    case 0
        mmObj.Model=[];
        mmObj.Order=0;
        mmObj=class(mmObj,'modelcomposite');
    case 1
        if iscell(varargin{1})
            for ii=1:length(varargin{1})
                mmObj.Model{ii}=varargin{1}{ii};
            end
            mmObj.Order=length(varargin{1});
        elseif isocmodel(varargin{1}) || ismultistagemodel(varargin{1})
            mmObj.Model{1}=varargin{1};
            mmObj.Order=1;
        else
        end
        mmObj=class(mmObj,'modelcomposite');

    otherwise
        counter=0;
        for ii=1:nargin
            if iscell(varargin{ii})
                for jj=1:length(varargin{ii})
                    counter=counter+1;
                    mmObj.Model{counter}=varargin{ii}{jj};
                end
            else
                counter=counter+1;
                mmObj.Model{counter}=varargin{ii};
            end
        end
        mmObj.Order=counter;
        mmObj=class(mmObj,'modelcomposite');
end